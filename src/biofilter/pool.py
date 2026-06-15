"""Parallel-pool abstraction shared by every heavy stage.

Locally we parallelise across CPU cores with ``multiprocessing``; on Sapelo2 we
scale across nodes with MPI via ``schwimmbad.MPIPool`` (which needs ``mpi4py``,
provided by the cluster's MPI module). The same driver code runs in both places:

    with get_pool(backend) as pool:
        sampler = emcee.EnsembleSampler(..., pool=pool)

MPIPool follows the controller/worker pattern -- worker ranks block in
``pool.wait()`` and exit, so callers must guard the science with
``if pool.is_master()`` when using MPI. :func:`get_pool` returns a uniform object
with ``.map``, ``is_master()`` and context-manager support for all backends.
"""
from __future__ import annotations

import contextlib
import multiprocessing as mp
import os


class _SerialPool:
    def map(self, fn, iterable):
        return list(map(fn, iterable))

    @staticmethod
    def is_master():
        return True

    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class _MPMaster:
    """multiprocessing wrapper exposing the same surface as MPIPool."""

    def __init__(self, nproc):
        # Force the 'fork' context: Python 3.14 defaults to 'forkserver' on
        # Linux, under which module-global state is not inherited. 'fork' keeps
        # the established behaviour; targets are also picklable classes, so this
        # is belt-and-suspenders.
        ctx = mp.get_context("fork")
        self._pool = ctx.Pool(nproc)

    def map(self, fn, iterable):
        return self._pool.map(fn, iterable)

    @staticmethod
    def is_master():
        return True

    def close(self):
        self._pool.close()
        self._pool.join()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


@contextlib.contextmanager
def get_pool(backend: str = "auto", nproc: int | None = None):
    """Yield a pool. ``backend`` in {auto, serial, multiprocessing, mpi}.

    ``auto`` picks MPI when launched under an MPI runtime (``OMPI_COMM_WORLD_SIZE``
    or ``PMI_SIZE`` set), otherwise multiprocessing.
    """
    if backend == "auto":
        backend = "mpi" if ("OMPI_COMM_WORLD_SIZE" in os.environ or "PMI_SIZE" in os.environ) else "multiprocessing"

    if backend == "serial":
        with _SerialPool() as p:
            yield p
        return

    if backend == "mpi":
        from schwimmbad import MPIPool  # imported lazily; needs mpi4py

        pool = MPIPool()
        if not pool.is_master():
            pool.wait()
            pool.close()
            raise SystemExit(0)        # worker ranks exit cleanly
        try:
            yield pool
        finally:
            pool.close()
        return

    # default: multiprocessing
    nproc = nproc or os.cpu_count()
    with _MPMaster(nproc) as p:
        yield p

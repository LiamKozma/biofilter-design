#!/bin/bash
# Biofilter-design Sapelo2 health diagnostics.
# Read-only: inspects queue, accounting, logs and outputs. Does NOT touch jobs.
# Usage (on a Sapelo2 login node):
#   bash slurm/diag.sh > diag_report4.txt 2>&1
ROOT="${BIOFILTER_ROOT:-/scratch/lmk04992/biofilter-design}"
RUN=sapelo2_full
N_DATASETS=400
cd "$ROOT" || { echo "FATAL: cannot cd to $ROOT"; exit 1; }

hr(){ printf '\n========== %s ==========\n' "$1"; }
# fatal patterns vs the benign scipy/numpy RuntimeWarnings we expect
FATAL='Traceback|Error:|error:|CANCELLED|DUE TO TIME LIMIT|oom|OOM|Out Of Memory|Killed|mutually exclusive|cpu bind|cpu-bind|ModuleNotFoundError|ImportError|No such file|Permission denied|srun: error|srun: fatal'
BENIGN='RuntimeWarning|invalid value encountered|divide by zero|overflow encountered|np.vstack|scipy/integrate|tag_match|UCX  WARN|Lmod|reloaded with a version'

echo "biofilter-design diagnostics  (diag_report4)"
echo "generated : $(date)"
echo "host      : $(hostname)"
echo "user      : $(whoami)"
echo "dir       : $(pwd)"

hr "1. QUEUE -- all my jobs"
squeue --me -o "%.20i %.18j %.2t %.11M %.11l %.4D %R"
hr "1b. my jobs grouped by state"
squeue --me -h -o "%t" | sort | uniq -c

# auto-detect the live job ids by name (newest of each)
SBC_ARR=$(squeue --me -h -n biofilt-sbc   -o "%A" | sort -un | tail -1)
CALIB=$(squeue --me   -h -n biofilt-calib -o "%A" | sort -un | tail -1)
SOBOL=$(squeue --me   -h -n biofilt-sobol -o "%A" | sort -un | tail -1)
echo
echo "detected job ids -> SBC array=${SBC_ARR:-none}  calib=${CALIB:-none}  sobol=${SOBOL:-none}"

hr "2. SBC ARRAY accounting (sacct -X, this array)"
if [ -n "$SBC_ARR" ]; then
  sacct -j "$SBC_ARR" -X --format=JobID%17,State%14,Elapsed,ExitCode,NodeList%10 2>/dev/null
  echo
  echo "--- task state counts ---"
  sacct -j "$SBC_ARR" -X -n --format=State 2>/dev/null | awk '{print $1}' | sort | uniq -c
  echo
  echo "--- REGRESSION CHECK: any bad terminal states? (want: none) ---"
  sacct -j "$SBC_ARR" -X -n --format=JobID%17,State,Elapsed 2>/dev/null \
     | grep -E 'TIMEOUT|FAILED|OUT_OF_ME|CANCELLED|NODE_FAIL' || echo "  none -- good"
  echo
  echo "--- peak memory of finished tasks (MaxRSS, top 10; limit = 24 x 2G = 48G/task) ---"
  sacct -j "$SBC_ARR" -n --format=JobID%20,MaxRSS,State 2>/dev/null | awk '$2 ~ /[0-9]/' | sort -k2 -h | tail -10
  echo
  echo "--- longest elapsed among finished tasks (limit 12:00:00) ---"
  sacct -j "$SBC_ARR" -X -n --format=Elapsed,State 2>/dev/null | grep -E 'COMPLETED|TIMEOUT|FAILED' | sort | tail -5
else
  echo "no SBC array currently in queue"
fi

hr "3. SBC OUTPUTS written so far"
nfiles=$(ls results/sbc_${RUN}_task*.json 2>/dev/null | wc -l)
echo "result files: ${nfiles} / ${N_DATASETS} datasets"
echo "--- 5 newest ---"
ls -t results/sbc_${RUN}_task*.json 2>/dev/null | head -5 | while read f; do
  printf '  %s  (%s bytes)\n' "$f" "$(stat -c%s "$f" 2>/dev/null)"
done

hr "4. Is a RUNNING SBC task actually PROGRESSING?"
if [ -n "$SBC_ARR" ]; then
  latest_out=$(ls -t logs/sbc_${SBC_ARR}_*.out 2>/dev/null | head -1)
  if [ -n "$latest_out" ]; then
    echo "most-recently-written .out: $latest_out  (mtime $(stat -c%y "$latest_out" 2>/dev/null))"
    n=$(grep -c 'dataset' "$latest_out" 2>/dev/null)
    echo "progress lines so far in this task: $n   (each = one finished dataset within the task)"
    echo "--- last progress lines ---"; grep 'dataset' "$latest_out" 2>/dev/null | tail -5
    echo "--- tail of .out ---"; tail -12 "$latest_out"
    latest_err="${latest_out%.out}.err"
    echo "--- .err: real errors only (benign RuntimeWarnings hidden) ---"
    grep -nE "$FATAL" "$latest_err" 2>/dev/null | grep -vE "$BENIGN" | tail -15 || true
    echo "    (benign RuntimeWarning lines suppressed: $(grep -cE "$BENIGN" "$latest_err" 2>/dev/null))"
  else
    echo "no sbc .out logs for array $SBC_ARR yet"
  fi
fi

hr "5. CALIBRATION job ($CALIB)"
if [ -n "$CALIB" ]; then
  echo "--- scontrol timing ---"
  scontrol show job "$CALIB" 2>/dev/null | grep -E 'RunTime|TimeLimit|JobState|NodeList' | sed 's/^ */  /'
  echo "--- calib .out tail (emcee progress / final 'acc=.. ESS~..' line) ---"
  tail -8 "logs/calib_${CALIB}.out" 2>/dev/null || echo "  no .out"
  echo "--- calib .err: real errors only ---"
  grep -nE "$FATAL" "logs/calib_${CALIB}.err" 2>/dev/null | grep -vE "$BENIGN" | tail -10 || true
  echo "--- calib .err last 3 raw lines (tqdm progress bar lives here) ---"
  tail -3 "logs/calib_${CALIB}.err" 2>/dev/null
fi
echo "--- calibration outputs (appear only when the job finishes) ---"
ls -la results/calib_chain_${RUN}.npy results/calib_summary_${RUN}.json 2>/dev/null || echo "  not written yet (job still running)"

hr "6. SOBOL output"
if [ -f results/sobol_${RUN}.json ]; then
  ls -la results/sobol_${RUN}.json
  ./.venv/bin/python -c "import json; d=json.load(open('results/sobol_${RUN}.json')); print('  top-level keys:', list(d)[:12])" 2>/dev/null \
    || echo "  (could not parse with venv python)"
else
  echo "  results/sobol_${RUN}.json not present"
fi

hr "7. RECENT log error scan (.err modified in last 3h)"
found=0
for f in $(find logs -name '*.err' -mmin -180 2>/dev/null); do
  real=$(grep -E "$FATAL" "$f" 2>/dev/null | grep -vcE "$BENIGN")
  if [ "${real:-0}" -gt 0 ]; then
    echo "  $f : $real real-error line(s)"
    grep -nE "$FATAL" "$f" 2>/dev/null | grep -vE "$BENIGN" | tail -3 | sed 's/^/      /'
    found=1
  fi
done
[ "$found" -eq 0 ] && echo "  no fatal errors in recently-modified logs -- good"

hr "8. DISK usage"
echo "--- sizes ---"
du -sh results logs 2>/dev/null
echo "--- scratch filesystem ---"
df -h . 2>/dev/null | tail -2
command -v quota >/dev/null 2>&1 && { echo "--- quota ---"; quota -s 2>/dev/null | tail -5; }

hr "9. ENVIRONMENT sanity (.venv)"
./.venv/bin/python -c "import sys,numpy,scipy,emcee; print('  python',sys.version.split()[0],'| numpy',numpy.__version__,'| scipy',scipy.__version__,'| emcee',emcee.__version__)" 2>&1 | head -3

hr "10. CONFIG echo (configs/${RUN}.yaml)"
grep -E 'run_name|n_walkers|n_steps|n_burn|n_base|n_datasets|inner_|backend|seed' configs/${RUN}.yaml 2>/dev/null | sed 's/^/  /'
echo "  sbc_array.sub --array line: $(grep -- '--array' slurm/sbc_array.sub 2>/dev/null | sed 's/#SBATCH//')"
echo "  sbc_array.sub --time  line: $(grep -- '--time'  slurm/sbc_array.sub 2>/dev/null | sed 's/#SBATCH//')"

hr "SUMMARY (heuristic)"
done_n=$(sacct -j "$SBC_ARR" -X -n --format=State 2>/dev/null | grep -c COMPLETED)
bad_n=$(sacct -j "$SBC_ARR" -X -n --format=State 2>/dev/null | grep -cE 'TIMEOUT|FAILED|OUT_OF_ME|NODE_FAIL')
run_n=$(squeue --me -h -t R -o '%t' | wc -l)
pend_n=$(squeue --me -h -t PD -o '%t' | wc -l)
echo "  SBC tasks COMPLETED        : ${done_n}"
echo "  SBC tasks TIMEOUT/FAILED   : ${bad_n}   $([ "${bad_n:-0}" -gt 0 ] && echo '<-- INVESTIGATE' || echo 'ok')"
echo "  SBC result files on disk   : ${nfiles} / ${N_DATASETS}"
echo "  jobs RUNNING / PENDING     : ${run_n} / ${pend_n}"
echo "  calib outputs present      : $([ -f results/calib_chain_${RUN}.npy ] && echo yes || echo 'not yet (running)')"
echo
echo "diagnostics complete -- $(date)"

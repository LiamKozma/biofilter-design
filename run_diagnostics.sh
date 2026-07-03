#!/bin/bash
cd /scratch/$USER/biofilter-design 2>/dev/null || cd ~/biofilter-design
echo "=== JOBS ==="; squeue --me -o "%.18i %.20j %.2t %.10M %.4D %R"
echo "=== which job is 46244325 ==="; scontrol show job 46244325 | grep -E "JobName|Command|RunTime|NodeList"
echo "=== recent logs ==="; ls -lt logs/ | head -30
echo "=== CALIB .err tail ==="; tail -n 25 logs/calib_*.err 2>/dev/null
echo "=== CALIB .out tail ==="; tail -n 15 logs/calib_*.out 2>/dev/null
echo "=== SOBOL .err tail ==="; tail -n 25 logs/sobol_*.err 2>/dev/null
echo "=== SOBOL .out tail ==="; tail -n 15 logs/sobol_*.out 2>/dev/null
echo "=== one SBC task (.out/.err) ==="; tail -n 15 logs/sbc_*_39.out 2>/dev/null; tail -n 15 logs/sbc_*_39.err 2>/dev/null
echo "=== results written so far ==="; ls -lt results/ | head; echo "sbc task files:"; ls results/sbc_sapelo2_full_task*.json 2>/dev/null | wc -l

  cd /scratch/$USER/biofilter-design
  echo "completed SBC task files:"; ls results/sbc_sapelo2_full_task*.json 2>/dev/null | wc -l
  echo "newest task .out (should show 'dataset N: ranks=...'):"
  tail -n 3 $(ls -t logs/sbc_*.out | head -1)
  squeue --me -h -t running -o "%i %M" | grep 46243364 | sort -k2 | tail -3


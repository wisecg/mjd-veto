#!/bin/bash
# This is a handy little tool for searching multiple patterns in the production
# log files, and figuring out the run they come from.
# C. Wiseman, 11/3/2016

# export logdir="/project/projectdirs/majorana/data/production/mjdProcessingLogs"
export logdir="./logs/DS0/"

# Event-level error checks ('s' denotes setting skip=true)
# s 1. Missing channels (< 32 veto datas in event)
# s 2. Extra Channels (> 32 veto datas in event)
# s 3. Scaler only (no QDC data)
#   4. Bad Timestamp: FFFF FFFF FFFF FFFF
# s 5. QDCIndex - ScalerIndex != 1 or 2
# s 6. Duplicate channels (channel shows up multiple times)
#   7. HW Count Mismatch (SEC - QEC != 1 or 2)
#   8. MJTRun run number doesn't match input file
# s 9. MJTVetoData cast failed (missing QDC data)
#   10. Scaler EventCount doesn't match ROOT entry
#   11. Scaler EventCount doesn't match QDC1 EventCount
#   12. QDC1 EventCount doesn't match QDC2 EventCount
# s 13. Indexes of QDC1 and Scaler differ by more than 2
# s 14. Indexes of QDC2 and Scaler differ by more than 2
#   15. Indexes of either QDC1 or QDC2 PRECEDE the scaler index
#   16. Indexes of either QDC1 or QDC2 EQUAL the scaler index
#   17. Unknown Card is present.
# s 18. Scaler & SBC Timestamp Desynch.
# s 19. Scaler Event Count reset.
# s 20. Scaler Event Count increment by > +1.
# s 21. QDC1 Event Count reset.
# s 22. QDC1 Event Count increment by > +1.
# s 23. QDC2 Event Count reset.
# s 24. QDC2 Event Count increment > +1.
# s 25. Buffer flush error.
#
# Run-level error checks
# 26. LED frequency very low/high, corrupted, or LED's off.
# 27. QDC threshold not found
# 28. No events above QDC threshold
# 29. Avg Panel LEDQDC deviates from expected mean by > 3 sigma.
# 30. nonLED Panel Hit Rate deviates from expected mean by 3 > sigma.

echo "For logs in $logdir -- "
var=$(grep -rI "Segmentation" $logdir | wc -l) && echo "Auto-veto segfaults: $var runs"
var=$(grep -rI "Sync failed" $logdir | wc -l) && echo "Sync issues: $var runs"
grep -rI "Sync failed" $logdir
var=$(grep -rI "Warning: missing start or stop" $logdir | wc -l) && echo "Duration issues: $var runs"
var=$(grep -rI "LED may be" $logdir | wc -l) && echo "LED issues: $var runs"
var=$(grep -rI "Error\[1\]:" $logdir | wc -l) && echo "E1 Missing Channels: $var runs"
var=$(grep -rI "Error\[4\]:" $logdir | wc -l) && echo "E4 Bad Timestamp: $var runs"
var=$(grep -rI "Error\[5\]:" $logdir | wc -l) && echo "E5 Index Mismatch: $var runs"
var=$(grep -rI "Error\[6\]:" $logdir | wc -l) && echo "E6 Duplicate Channels: $var runs"
var=$(grep -rI "Error\[13\]:" $logdir | wc -l) && echo "E13 QDC1 Index Mismatch: $var runs"
var=$(grep -rI "Error\[14\]:" $logdir | wc -l) && echo "E14 QDC2 Index Mismatch: $var runs"
var=$(grep -rI "Error\[18\]:" $logdir | wc -l) && echo "E18 Scaler/SBC Desynch: $var runs"
var=$(grep -rI "Error\[19\]:" $logdir | wc -l) && echo "E19 SEC Reset: $var runs"
var=$(grep -rI "Error\[20\]:" $logdir | wc -l) && echo "E20 SEC Jump: $var runs"
var=$(grep -rI "Error\[21\]:" $logdir | wc -l) && echo "E21 QDC1 EC reset: $var runs"
var=$(grep -rI "Error\[22\]:" $logdir | wc -l) && echo "E22 QDC1 EC jump: $var runs"
var=$(grep -rI "Error\[23\]:" $logdir | wc -l) && echo "E23 QDC2 EC reset: $var runs"
var=$(grep -rI "Error\[24\]:" $logdir | wc -l) && echo "E24 QDC2 EC jump: $var runs"
var=$(grep -rI "Error\[25\]:" $logdir | wc -l) && echo "E25 Buffer Flush: $var runs"
var=$(grep -rI "Error\[29\]:" $logdir | wc -l) && echo "E29 LED QDC Deviation: $var runs"
var=$(grep -rI "Error\[30\]:" $logdir | wc -l) && echo "E30 Hit Rate Deviation: $var runs"

# just show the output
# grep -rI "Error\[18\]:" $logdir

# Output two lines from the same file with an OR:
# echo "LED issues " && grep -rIl "LED may be" $logdir | xargs grep -hn "processing run\|LED may be" # | wc -l

# how many hits?
# hitType="2+ panels"
# hitType="vertical"
# hitType="side+bottom"
# hitType="top+sides"
# hitType="compound"
# grep -rIl "$hitType" $logdir | xargs grep -hn "$hitType" | wc -l



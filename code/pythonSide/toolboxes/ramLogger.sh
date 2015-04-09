#!/bin/bash

clear

echo "Starting RAM Logging now"
echo 
echo "output log full path provided = $1"
echo "sleep time will be $2 seconds"

if [ -f $1 ]; then
	rm $1 
	echo "removed previous file: $1"
fi

echo "Start of log file:" >>$1
echo "sleep time will be $2 seconds">>$1
echo "-----------------------------">>$1
echo "">>$1


while true; do 
			if [ -f $1 ]; then
				free -m >> $1
				echo "Just logged latest usage, after pressing Ctrl+c to clean and plot type:"
				echo "'python /mnt/Data1/Todai_Work/EclipseWorkspace/SMODT/code/pythonSide/toolboxes/memLogCleanAndPlot.py'"
				sleep $2
			fi

			if [ ! -f $1 ]; then
				break
			fi

			done

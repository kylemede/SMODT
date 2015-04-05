#!/bin/bash

clear

echo "Starting RAM Logging now"
echo 
echo "output log full path provided = $1"
echo "sleep time will be $2 seconds"

rm $1

while true; do 
			free -m >> $1
			echo "Just logged latest usage, after pressing Ctrl+c to clean and plot type:"
			echo "'python /mnt/Data1/Todai_Work/EclipseWorkspace/SMODT/code/pythonSide/toolboxes/memLogCleanAndPlot.py'"
			sleep $2
			done

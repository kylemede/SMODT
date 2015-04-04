#!/bin/bash

clear

echo "Starting RAM Logging now"
echo 

rm /mnt/Data1/Todai_Work/Data/data_SMODT/RAMusage.log

while true; do 
			free -m >> /mnt/Data1/Todai_Work/Data/data_SMODT/RAMusage.log
			echo "Just logged latest usage, after pressing Ctrl+c to clean and plot type:"
			echo "'python /mnt/Data1/Todai_Work/EclipseWorkspace/SMODT/code/pythonSide/toolboxes/memLogCleanAndPlot.py'"
			sleep 600
			done

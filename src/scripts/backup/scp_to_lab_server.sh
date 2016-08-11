#!/bin/bash
# Shell script to copy CollecTF backup to the lab server.

ERRORFILE="/home/erilllab/backups/backup.err"

MOST_RECENT_SQL=$(ls -Art /home/erilllab/backups/*.sql.gz | tail -n -1)

scp $MOST_RECENT_SQL erilllab@erilllab.biosci.umbc.edu:./collectf_backups/collectf_backup.sql.gz 2>>$ERRORFILE

if [ $? -eq 0 ]
then
    msg="[$(date)] mysqldump copied succesfully."
else
    msg="[$(date)] mysqldump copied error. See $ERRORFILE."
fi

# send the email
echo $msg | mail -s "CollecTF: $msg" sefa1@umbc.edu

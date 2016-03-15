#!/bin/bash
# Shell script to backup CollecTF MySQL database.

MyUSER=""                   # USERNAME
MyPASS=""               # PASSWORD
MyHOST="localhost"              # Hostname
DB="collectfdb"

# Linux bin paths
MYSQL="$(which mysql)"
MYSQLDUMP="$(which mysqldump)"
CHOWN="$(which chown)"
CHMOD="$(which chmod)"
GZIP="$(which gzip)"
FIND="$(which find)"

# Backup Dest directory, change this if you have someother location
DEST="/home/erilllab/backups"

# Main directory where backup will be stored
MBD="$DEST/"

# Get data in dd-mm-yyyy format
NOW="$(date +"%m-%d-%Y-%H-%M-%S")"

# Backup here
FILE="$MBD/collectfdb.$NOW.sql.gz"
LOGFILE="$MBD/backup.log"
ERRORFILE="$MBD/backup.err"

$MYSQLDUMP $DB 2>>$ERRORFILE | $GZIP -1 >$FILE

if [ ${PIPESTATUS[0]} -eq 0 ]
then
    echo "[$NOW] mysqldump completed succesfully." >>$LOGFILE
    tac $LOGFILE | mail -s "[$NOW] collectf mysqldump success." sefa1@umbc.edu
else
    echo "[$NOW] mysqldump encountered a problem. Look in $ERRORFILE for information" >>$LOGFILE
    tac $LOGFILE | mail -s "[$NOW] collectf mysqldump error" sefa1@umbc.edu
fi

# remove files older than 7 days
$FIND $MBD/*.sql.gz -type f -mtime +7 | xargs rm -f

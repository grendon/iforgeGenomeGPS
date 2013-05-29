#!/bin/sh

echo "----------------------------------------"
echo "Epilogue Args:"
echo "Job ID: $1"
echo "User ID: $2"
echo "Group ID: $3"
echo "Job Name: $4"
echo "Session ID: $5"
echo "Resource List: $6"
echo "Resources Used: $7"
echo "Queue Name: $8"
echo "Account String: $9"
echo "Exit status: ${10}"
echo "----------------------------------------"

echo "----------------------------------------"
echo "checkjob output:"
/usr/local/moab/bin/checkjob $1
echo "----------------------------------------"

exit 0

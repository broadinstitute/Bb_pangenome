#!/bin/bash

# Michael J. Foster
# github.com/mjfos2r
# pre-commit hook to automatically track large files with LFS

# Max filesize: 100MB (per github)
MAX_SIZE_MB=100
MAX_SIZE_BYTES=$((MAX_SIZE_MB * 1024 * 1024))

echo "Checking files to ensure that all are below ${MAX_SIZE_MB}MB... please stand by..."

# Get all of the files we've staged in this commit
STAGED_FILES=$(git diff --cached --name-only)
LFS_MODIFIED=0

for FILE in $STAGED_FILES; do
    if [ ! -f "$FILE" ]; then
		    continue
		fi

    FILE_SIZE=$(stat -f%z "$FILE" 2>/dev/null || stat -c%s "$FILE" 2>/dev/null)
		FILE_SIZE=$((FILE_SIZE + 0)) # force to integer?
    FILE_SIZE_MB=$(echo "scale=2; $FILE_SIZE / 1048576" | bc) # format in MB for human readability.

    # now check that it's large enough to use LFS.
    if [ "$FILE_SIZE" -gt "$MAX_SIZE_BYTES" ]; then
			  echo "Large file detected: $FILE (${FILE_SIZE_MB}MB)"
				FILE_EXTENSION="${FILE##*.}"

				# check if it's already tracked by LFS
				if ! grep -a "\"*.$FILE_EXTENSION\"" .gitattributes 2>/dev/null; then
					echo "Setting up LFS tracking for *.$FILE_EXTENSION files"
					git lfs track "*.$FILE_EXTENSION"
					git add .gitattributes
					LFS_MODIFIED=1
				fi
				# if file is staged but not as LFS, re-stage it.
				IS_LFS=$(git check-attr filter "$FILE" | grep -q "lfs" && echo "true" || echo "false")
				if [ "IS_LFS" = "false" ]; then
					echo "Re-staging $FILE to use LFS"
					git reset -q HEAD "$FILE"
					git add "$FILE"
					LFS_MODIFIED=1
				fi
    fi
done

if [ $LFS_MODIFIED -eq 1  ]; then
	echo "Large files have been configured to use Git LFS"
	echo "Make sure LFS is set up in this repository."
fi

exit 0



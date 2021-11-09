#!usr/bin/bash

# Records have the form "location/name/date"
# We're only going to search for "name/date", without location
awk -F '/' '{ print $2"/"$3}' pangomiss.txt > pangomiss_noloc.txt

# If file exists, re-initialize it
touch pangoalts.txt && rm pangoalts.txt && touch pangoalts.txt


pangomiss="pangomiss_noloc.txt"
lines=$(cat $pangomiss)

for line in $lines
do
    echo $line
    thisgrep=`xzgrep $line sequences_metadata.txt`
    # If there are more than one match, add to the same line
    var=""
    for element in $thisgrep
    do
        var+="${element}"
        var+="," # Not the greatest formatting, but will work
    done
    echo "$line, $var" >> pangoalts.txt
done

head -20 pangoalts.txt


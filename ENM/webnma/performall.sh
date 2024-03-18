#!/bin/bash

# Define an array of NMR files
# NMRfiles=('1A03' '2LGV' '7UO6' '7R0R' '2KB1' '7VIL')
NMRfiles=('8F4V' '7ZE0' '4TRX' '1CIS' '2ABD')
# Loop through each NMR file in the array
for file in "${NMRfiles[@]}"
do
  # Perform a function on each NMR file
  # In this example, we'll just print the filename to the console
  echo "Processing NMR file $file"
  if [ ! -d $file ]; then
  # Create the directory if it doesn't exist
  	mkdir -p $file
  fi
  webnma ca -s -p "$file" "$file"_*
done

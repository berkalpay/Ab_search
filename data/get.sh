set -e

subject_re="supplement/(.*)_H"
tarfn_re="supplement/(.*)"
while read u; do
    # Read meta-info
    if [[ $u =~ $subject_re ]]
    then
	subject="${BASH_REMATCH[1]}"
    fi
    echo $subject
    if [[ $u =~ $tarfn_re ]]
    then
	tarfn="${BASH_REMATCH[1]}"
    fi

    # Download
    wget -nv $u

    # Unzip
    tar -xvzf "${tarfn/\%2B/+}"
    rm "${tarfn/\%2B/+}"

    # Remove unwanted columns
    mkdir $subject
    for fn in consensus-cdr3nt*_minimal/*.txt; do
	# Extract replicate name
	rep_re="minimal/(.*)_consensus.txt"
	if [[ $fn =~ $rep_re ]]
	then
	    rep="${BASH_REMATCH[1]}"
	fi

	# Create trimmed file
	awk -F',' '{OFS=","; print $3, $4, $5, $7, $9, $13, $25}' $fn > "${subject}/${rep}.csv"
    done

    # Delete raw files
    rm -r consensus-cdr3nt*_minimal
done < urls.txt

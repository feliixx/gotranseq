#!/bin/bash
# Generate data for tests 
# to use it, run 
#   ./generateTestData.sh test.fna 


# array of parameters to test
tests=("-frame=1"
"-frame=2"
"-frame=3"
"-frame=F"
"-frame=-1"
"-frame=-2"
"-frame=-3"
"-frame=R"
"-frame=6"
"-frame=6 -clean"
"-frame=6 -trim"
"-frame=6 -alternative"
"-frame=6 -table=2"
"-frame=6 -table=3"
"-frame=6 -table=4"
"-frame=6 -table=5"
"-frame=6 -table=6"
"-frame=6 -table=9"
"-frame=6 -table=10"
"-frame=6 -table=11"
"-frame=6 -table=12"
"-frame=6 -table=13"
"-frame=6 -table=14"
"-frame=6 -table=16"
"-frame=6 -table=21"
"-frame=6 -table=22")



echo '[' > data.json

for test in "${tests[@]}"
do 
    echo "{" >> data.json
    echo "\"options\": \"$test\"," >> data.json
    # generate the expected output
    transeq -sequence $1 -outseq out.faa ${test[0]} ${test[1]}
    echo -n "\"expected\": \""  >> data.json 
    while IFS= read -r line 
      do 
        echo -n $line"\n" >> data.json
    done < out.faa
    echo  "\"" >> data.json
    echo "}," >> data.json

done

# remove the last ',' to get a valid json
sed -i '$ s/.$//' data.json
echo ']' >> data.json
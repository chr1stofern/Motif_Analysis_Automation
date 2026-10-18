# Compare matching .tsv files between dir1 and dir2
for file in ./data/*.tsv; do
    filename=$(basename "$file")
    if [ -f "intermediate_files/preprocessed_sequences/$filename" ]; then
        echo "=== Checking differences for $filename ==="
        diff -u "$file" "intermediate_files/preprocessed_sequences/$filename"
    fi
done

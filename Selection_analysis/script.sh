#For orthologues downloaded from NCBI/Ensembl
base=MYD88 #Gene name

#RefSeq transcript edit
sed 's/^>.*\[/>/g' -i ${base}_refseq_protein.fasta
sed 's/\]//g' -i ${base}_refseq_protein.fasta
sed 's/ /_/g' -i ${base}_refseq_protein.fasta
sed 's/ /_/g' -i ${base}_refseq_protein.fasta
cat ${base}_refseq_protein.fasta | cut -d"_" -f1,2 > tmp
mv tmp ${base}_refseq_protein.fasta

###check the headers
grep ">" ${base}_refseq_protein.fasta

cat  ${base}_refseq_transcript.fasta | sed 's/^>.*: />/g' | cut -d" " -f1,2 >  ${base}_refseq_transcript.fna
sed 's/ /_/g' -i ${base}_refseq_transcript.fna
sed '/^[[:space:]]*$/d' -i ${base}_refseq_transcript.fna
cat *.fna | cut -d"_" -f1,2 > tmp
 mv tmp ${base}_refseq_transcript.fna
#Check the fasta header
grep ">" ${base}_refseq_transcript.fna 


pal2nal.pl ${base}_refseq_protein.aln.fasta ${base}_refseq_transcript.fna -output fasta -nogap > ${base}_refseq_codon.aln.fasta

#Phylogenetic tree
clustalo -i  ${base}_refseq_protein.fasta -o ${base}_refseq_protein.aln.fasta
iqtree2 -s ${base}_refseq_protein.aln.fasta -T 20 -v

#Edit tree to include {fg} format of branch selection
cp ${base}_refseq_protein.aln.fasta.treefile ${base}_refseq_protein.aln.fasta.treefile_fg
vi ${base}_refseq_protein.aln.fasta.treefile #{fgg}, {fga}, {fgco}, {fgca}
vi ${base}_refseq_protein.aln.fasta.treefile_fg #{fg} for all four species

hyphy BUSTED --alignment ${base}_refseq_codon.aln.fasta --tree ${base}_refseq_protein.aln.fasta.treefile --output ${base}.busted.json
hyphy BUSTED --alignment ${base}_refseq_codon.aln.fasta --tree ${base}_refseq_protein.aln.fasta.treefile_fg --output ${base}.busted.json --branches fg

for i in fga fgco fgca fgg
do
for p in FEL MEME 
do
hyphy ${p} --alignment ${base}_refseq_codon.aln.fasta --tree ${base}_refseq_protein.aln.fasta.treefile --branches $i --output ${i}.${p}.json
done
done

mpiexec -n 20 HYPHYMPI CPU=20 FUBAR --alignment TLR7_refseq_codon.aln.fasta.best-gard --branches fg --output ${base}.FUBAR.json

seqkit grep -p "Cygnus_atratus" ${base}_refseq_protein.fasta
seqkit grep -p "Anas_platyrhynchos" ${base}_refseq_protein.fasta
seqkit grep -p "Cygnus_olor" ${base}_refseq_protein.fasta
seqkit grep -p "Gallus_gallus" ${base}_refseq_protein.fasta

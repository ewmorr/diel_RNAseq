
grep "K12313" P1P2P3.KO | cut -f 1 > K12313_gene_ids.txt

perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' K12313_gene_ids.txt genes_P1P2P3.derep.fna > K12313_genes.fa

grep "K10355" P1P2P3.KO | cut -f 1 > K10355_gene_ids.txt

perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' K10355_gene_ids.txt genes_P1P2P3.derep.fna > K10355_genes.fa

grep "K12314" P1P2P3.KO | cut -f 1 > K12314_gene_ids.txt

perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' K12314_gene_ids.txt genes_P1P2P3.derep.fna > K12314_genes.fa


metagenematrix = load('sample_data/example_metagene3.txt');

TF = ThreeFactorsMutualInhibitionAndCofactor(500,100,1);

[GE,MetageneMatrix] = TF2GE(TF,'MetaGeneMatrix',metagenematrix);
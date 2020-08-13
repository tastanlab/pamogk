function [indx]= normalized_kmeans(U,numclass)
stream = RandStream.getGlobalStream;
reset(stream);
U_normalized = U ./ repmat(sqrt(sum(U.^2, 2)), 1,numclass);

indx = litekmeans(U_normalized,numclass, 'MaxIter',100, 'Replicates', 50);

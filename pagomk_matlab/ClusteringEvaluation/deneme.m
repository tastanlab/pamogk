indx = indx(:);
[newIndx] = bestMap(Y,indx);
res(1) = mean(Y==newIndx);
res(2) = MutualInfo(Y,newIndx);
res(3) = purFuc(Y,newIndx);
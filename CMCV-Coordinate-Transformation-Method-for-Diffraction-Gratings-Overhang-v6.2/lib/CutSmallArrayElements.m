function A=CutSmallArrayElements(A,accuracy)

idx = abs(A) < accuracy;
A(idx) = 0;

end
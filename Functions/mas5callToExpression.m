function [expression,expressionid]=mas5callToExpression(calls,dataset_probe,probes,probe2geneid,modelGenes)
%INPUT:
%calls: the output of preprocessing (mas5 call)
%dataset_probe:probe in the data set
%probes:all probes of the platform with corresponding gene id
%probe2geneid:corresponding gene ids of the probes. probes and probe2geneid
%come from the probes2id.mat
%modelGenes:genes in the generic metabolic model
%OUTPUT:
%expression:a column vector of how often a gene is called present across
%samples. Value between 0 and 1 (inclusive)
%expressionid:corresponding gene ids of the expression values.
[c,ind]=ismember(dataset_probe,probes);
expressionid=probe2geneid(ind(c),1);
calls=calls(c,:);
c=ismember(expressionid,modelGenes);
expression=calls(c,:);
expressionid=expressionid(c,:);
[expression,expressionid]=uniq_exp(expression,expressionid);
expression=expression>=2;
expression=sum(expression,2)/size(expression,2);
GB comments:
Prob1: 100%
Prob2:
P1:100 
P2: 50 the code does not consider whether the ORF is in correct reading frame. 
P3:100
P4:100
P5:100 
Prob3
P1: 100 
P2:100
P3:0 
Overall: 83


% Homework 1. Due before class on 9/5/17

%% Problem 1 - addition with strings

% Fill in the blank space in this section with code that will add 
% the two numbers regardless of variable type. Hint see the matlab
% functions ischar, isnumeric, and str2num. 

%your code should work no matter which of these lines is uncommented. 
x = 3; y = 5; % integers

x = '3'; y= '5'; %strings

%x = 3; y = '5'; %mixed

%your code goes here

if ischar(x)
    x = str2num(x)
end
if ischar(y)
    y = str2num(y)
end

z = x + y


%output your answer

%% Problem 2 - our first real biology problem. Open reading frames and nested loops.

%part 1: write a piece of code that creates a random DNA sequence of length
% N (i.e. consisting of the letters ATGC) where we will start with N=500 base pairs (b.p.).
% store the output in a variable
% called rand_seq. Hint: the function randi may be useful. 
% Even if you have access to the bioinformatics toolbox, 
% do not use the builtin function randseq for this part. 

N = 500; % define sequence length

sequence = "";
rand_seq = randi([1 4],[1,N]);

for i = 1:N
    if rand_seq(i) == 1
        sequence = sequence + "T";
    end
    if rand_seq(i) == 2
        sequence = sequence + "A";
    end
    if rand_seq(i) == 3
        sequence = sequence + "C";
    end
    if rand_seq(i) == 4
        sequence = sequence + "G";
    end
end
%part 2: open reading frames (ORFs) are pieces of DNA that can be
% transcribed and translated. They start with a start codon (ATG) and end with a
% stop codon (TAA, TGA, or TAG). Write a piece of code that finds the longest ORF 
% in your seqeunce rand_seq. Hint: see the function strfind.

startcodon = strfind(sequence, 'ATG');
stopcodon1 = strfind(sequence, 'TAA');
stopcodon2 = strfind(sequence, 'TGA');
stopcodon3 = strfind(sequence, 'TAG');

x = cat(2, stopcodon1,stopcodon2,stopcodon3);

for i = 1:length(startcodon)
    length_orf(i) = min(x(startcodon(i) < x)) - startcodon(i);
end

%part 3: copy your code in parts 1 and 2 but place it inside a loop that
% runs 1000 times. Use this to determine the probability
% that an sequence of length 500 has an ORF of greater than 50 b.p.

count = 0;
N = 500; % define sequence length
for k = 1:1000
sequence = "";
rand_seq = randi([1 4],[1,N]);
    for i = 1:N
        if rand_seq(i) == 1
            sequence = sequence + "T";
        end
        if rand_seq(i) == 2
            sequence = sequence + "A";
        end
        if rand_seq(i) == 3
            sequence = sequence + "C";
        end
        if rand_seq(i) == 4
            sequence = sequence + "G";
        end
    end
   
    startcodon = strfind(sequence, "ATG");
    stopcodon1 = strfind(sequence, "TAA");
    stopcodon2 = strfind(sequence, "TGA");
    stopcodon3 = strfind(sequence, "TAG");

    x = cat(2, stopcodon1,stopcodon2,stopcodon3);

    for i = 1:length(startcodon)
        if min(x(startcodon(i) < x)) ~= 0
            length_orf(i) = (min(x(startcodon(i) < x)) - startcodon(i));
            if max(length_orf) > 50
                count = count + 1;
            end
        end
    end
   
end
probability = (count/1000);


%part 4: copy your code from part 3 but put it inside yet another loop,
% this time over the sequence length N. Plot the probability of having an
% ORF > 50 b.p. as a funciton of the sequence length. 


sequence_length = 500; % define sequence length

for N = 1:length(sequence_length)
count = 0;
N = 500; % define sequence length
for k = 1:1000
sequence = "";
rand_seq = randi([1 4],[1,N]);
    for i = 1:N
        if rand_seq(i) == 1
            sequence = sequence + "T";
        end
        if rand_seq(i) == 2
            sequence = sequence + "A";
        end
        if rand_seq(i) == 3
            sequence = sequence + "C";
        end
        if rand_seq(i) == 4
            sequence = sequence + "G";
        end
    end
   
    startcodon = strfind(sequence, "ATG");
    stopcodon1 = strfind(sequence, "TAA");
    stopcodon2 = strfind(sequence, "TGA");
    stopcodon3 = strfind(sequence, "TAG");

    x = cat(2, stopcodon1,stopcodon2,stopcodon3);

    for i = 1:length(startcodon)
        if min(x(startcodon(i) < x)) ~= 0
            length_orf(i) = (min(x(startcodon(i) < x)) - startcodon(i));
            if max(length_orf) > 50
                count = count + 1;
            end
        end
    end
end
probability(N) = (count/1000);
end

figure;
plot(probability);


%part 5: Make sure your results from part 4 are sensible. What features
% must this curve have (hint: what should be the value when N is small or when
% N is very large? how should the curve change in between?) Make sure your
% plot looks like this. 

When the N is small, there shoul be a lower probability that there are ORFs that are greater than 50bps vs. when there is a higher N where that will increase.

%% problem 3 data input/output and simple analysis

%The file qPCRdata.txt is an actual file that comes from a Roche
%LightCycler qPCR machine. The important columns are the Cp which tells
%you the cycle of amplification and the position which tells you the well
%from the 96 well plate. Each column of the plate has a different gene and
%each row has a different condition. Each gene in done in triplicates so
%columns 1-3 are the same gene, columns 4-6 the same, etc.
%so A1-A3 are gene 1 condition 1, B1-B3 gene 1 condition 2, A4-A6 gene 2
%condition 1, B4-B6 gene2 condition 2 etc. 

% part1: write code to read the Cp data from this file into a vector. You can ignore the last two
% rows with positions beginning with G and H as there were no samples here. 

T = readtable('qPCRdata.txt');

T1 = T(1:72, 5:5);

T2 = table2array(T1);

% Part 2: transform this vector into an array representing the layout of
% the plate. e.g. a 6 row, 12 column array should that data(1,1) = Cp from
% A1, data(1,2) = Cp from A2, data(2,1) = Cp from B1 etc. 
% 
B = reshape (T2, [6,12])

% Part 3. The 4th gene in columns 10 - 12 is known as a normalization gene.
% That is, it's should not change between conditions and it is used to normalize 
% the expression values for the others. For the other three
% genes, compute their normalized expression in all  conditions, normalized to condition 1. 
% In other words, the fold change between these conditions and condition 1. The
% formula for this is 2^[Cp0 - CpX - (CpN0 - CpNX)] where Cp0 is the Cp for
% the gene in the 1st condition, CpX is the value of Cp in condition X and
% CpN0 and CpNX are the same quantitites for the normalization gene.
% Plot this data in an appropriate way. 



%% Challenge problems that extend the above (optional)

% 1. Write a solution to Problem 2 part 2 that doesn't use any loops at
% all. Hint: start by using the built in function bsxfun to make a matrix of all distances
% between start and stop codons. 

% 2. Problem 2, part 4. Use Matlab to compute the exact solution to this
% problem and compare your answer to what you got previously by testing
% many sequences. Plot both on the same set of axes. Hint: to get started 
% think about the following:
% A. How many sequences of length N are there?
% B. How many ways of making an ORF of length N_ORF are there?
% C. For each N_ORF how many ways of position this reading frame in a
% sequence of length N are there?

% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty



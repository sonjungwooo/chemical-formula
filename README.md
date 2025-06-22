# chemical_formula 
> Made for sci-competition
> | dedicate to YJ

## Research Topic

Developing an app using programming that automatically completes chemical equations to help students memorize them more easily and efficiently.

⸻

## Motivation

Students in our age group typically learn about molecules and atoms in the 8th grade and begin memorizing chemical reactions more seriously in the 9th grade. However, many of our classmates found it difficult to memorize these reactions, and we also experienced confusion when we first encountered them. To make the memorization process easier and more accessible, we decided to create a program that helps users understand and recall chemical equations more effectively.

⸻

## Preliminary Research and Considerations

A chemical equation represents the substances involved in a chemical reaction using symbols and formulas. When a reaction occurs, the reactants undergo atomic rearrangement, forming new bonds and resulting in products. One such bonding process is ionic bonding, which is mainly taught in middle school. However, in addition to ionic bonding, there are also metallic and covalent bonds.

When writing a chemical equation, the charges of the atoms involved are not always balanced. To address this, the number of atoms must be adjusted to satisfy the law of conservation of mass. This includes balancing the coefficients of each molecule in the equation. While middle school textbooks primarily focus on ionic bonding, our program aims to cover various types of chemical bonds and support automatic balancing of chemical equations.

⸻

## Procedure
	1.	We first considered which software and programming language to use.
	2.	We designed the program logic before coding:
	•	Accept a chemical reaction input from the user
	•	Parse and split the chemical formula
	•	Analyze the atoms involved and determine the correct coefficient ratios
	3.	Once the program structure was finalized, we began coding.
	4.	During development, we encountered and fixed various bugs through trial and error.
	5.	After initial debugging, we tested the program repeatedly to identify and correct any remaining issues.

⸻

## Methodology

We used the C++ programming language and developed the program using Visual Studio. By frequently running the program, we continuously checked for errors. To balance the coefficients of the chemical equations, we applied the method of undetermined coefficients.

Each coefficient was treated as an unknown variable (e.g., n1, n2, n3, n4), and we created a system of linear equations based on the number of atoms on both sides of the reaction. Using one variable as a reference (e.g., n1), we derived the relative values of the others and then simplified them into the smallest possible ratio of natural numbers.

To ensure the robustness of the program, we tested it with various chemical reactions. By inputting equations with different structures and complexities, we were able to find and resolve additional errors, gradually improving the program’s accuracy and reliability.

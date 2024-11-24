# Project Aim
The aim of this project is to use Mutation testing to test a real-world software project
with the help of open-source tools.

# What is Mutation Testing?
Mutation testing is a software testing technique that assesses the effectiveness of a test suite by introducing small, intentional changes (mutations) into the source code and determining whether the existing tests can detect these changes. Each mutation represents a potential bug, and the goal is to ensure that the test suite is sensitive enough to identify such alterations. If a mutation is not detected by the tests, it implies a weakness in the test suite, indicating that the code coverage or the quality of the tests may be insufficient. Mutation testing helps developers identify areas of improvement in their test suites and enhances overall software reliability by providing insights into the robustness of the testing process. <br/>

Mutation testing involves two paradigms for evaluating the effectiveness of test cases: Weak Mutation Killing and Strong Mutation Killing.
* In the Weak Mutation Killing paradigm, the memory state of the program after the execution of a mutated statement differs from the memory state when the statement remains unaltered and is executed. Interestingly, in this scenario, the program's output on a given test case may remain unchanged, regardless of whether a program statement has undergone mutation.
* Conversely, in the Strong Mutation Killing paradigm, the output of the program on a given test case must demonstrate a noticeable change when a statement undergoes mutation compared to its unaltered counterpart. Importantly, strong killing of mutants indicates that the error introduced by mutation propagates through the program, leading to distinct outputs in the presence and absence of the mutant.

Our adoption of mutation testing as a testing strategy is characterized by a deliberate emphasis on the robust neutralization of mutants through strong mutation killing.
This strategic approach aims to thoroughly validate the resilience and accuracy of the software under examination, ensuring that mutations result in discernible changes in the program's behavior and output.

# Testing Strategies and Tools used
Chosen testing strategy mainly focuses on strong mutation killing. This type of killing indicates that the error introduced by the mutation propagates through the program.
Tools used:
* VS Code IDE
* Pitest Java Library
* Maven
* JUnit

# Results
Total Lines in the code= 1440 (App.java) + 520(AppTest.java)=1960 <br/>
Our final test case design produced excellent results for us i.e., Line Coverage - 97% and Mutation Coverage - 81%. <br/>
Initially, the scores were a bit low, but by examining the mutants that survived and creating test cases to raise the mutation coverage percentage, we attempted to increase this mutation coverage. We examined the lines that the test cases did not cover and added test cases accordingly. <br/>

# Contributions
We had a long discussion regarding the project before deciding to perform mutation testing on the string-algorithms codebase. Given our mutual familiarity with Java, we chose to write the algorithm implementations in that language. After dividing the functions equally among ourselves, we eventually combined the code. We additionally provided test cases that would eliminate the mutants for each implementation of the function.
* Analysis and Testcases for Algo 1-14 - Arin Awasthi
* Analysis and Testcases for Algo 15-28 - Jainav Sanghvi



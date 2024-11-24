package st;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.junit.Test;

/**
 * Unit test for Array Functions Library.
 */
public class AppTest {
    App obj = new App();
    @Test
    public void TestRabinKarp() {
        String txt = "ABCFGHIJKLMNOPQRSTUVWXZXYZOPQRSTUWXYZ";
        String pat = "XYZOPQRS";
        int q = 101;
        // Check that the pattern "XYZOPQRS" starts at index 23
        assertEquals(23, obj.rabinKarp(pat, txt, q));

        txt = "ABCDXYEFGHIJK";
        pat = "ABCDX";
        // Pattern "ABCDX" is found at the beginning (index 0)
        assertEquals(0, obj.rabinKarp(pat, txt, q));

        txt = "KLMNOPQRSTUVWXZABCFGHIJKLMNOP";
        pat = "XYZOPQRSUVW";
        // Pattern "XYZOPQRSUVW" does not appear in the text, should return -1
        assertEquals(-1, obj.rabinKarp(pat, txt, q));

        txt = "ASDFGH";
        pat = "GHI";
        // Pattern "GHI" is not found in the text, should return -1
        assertEquals(-1, obj.rabinKarp(pat, txt, q));

        txt = "XZXZXZXZXZXZXZXY";
        pat = "XZXZXY";
        // Pattern "XZXZXY" starts at index 10
        assertEquals(10, obj.rabinKarp(pat, txt, q));
    }

    @Test
    public void TestZAlgorithm() {
        // Test Case 1: Pattern exists at the start of the text
        String txt = "ABCDXYEFGHIJK";
        String pat = "ABCDX";
        assertEquals(0, obj.ZAlgorithm(txt, pat));  // Pattern "ABCDX" starts at index 0 in the text.

        // Test Case 2: Pattern doesn't exist in the text
        txt = "KLMNOPQRSTUVWXZABCFGHIJKLMNOP";
        pat = "XYZOPQRSUVW";
        assertEquals(-1, obj.ZAlgorithm(txt, pat));  // Pattern doesn't exist in the text.

        // Test Case 3: Pattern is not found in the text (no match)
        txt = "ASDFGH";
        pat = "GHI";
        assertEquals(-1, obj.ZAlgorithm(txt, pat));  // Pattern "GHI" is not found in the text.

        // Test Case 4: Pattern occurs at the end of the text
        txt = "XZXZXZXZXZXZXZXY";
        pat = "XZXZXY";
        assertEquals(10, obj.ZAlgorithm(txt, pat));  // Pattern "XZXZXY" starts at index 10 in the text.
    }


    @Test
    public void TestKMPAlgorithm() {
        String txt = "ABCFGHIJKLMNOPQRSTUVWXZXYZOPQRSTUWXYZ";
        String pat = "XYZOPQRS";
        assertEquals(23, obj.KMPSearch(pat, txt), 0.0);

        txt = "ABCDXYEFGHIJK";
        pat = "ABCDX";
        assertEquals(0, obj.KMPSearch(pat, txt), 0.0);

        txt = "KLMNOPQRSTUVWXZABCFGHIJKLMNOP";
        pat = "XYZOPQRSUVW";
        assertEquals(-1, obj.KMPSearch(pat, txt), 0.0);

        txt = "ASDFGH";
        pat = "GHI";
        assertEquals(-1, obj.KMPSearch(pat, txt), 0.0);

        txt = "XZXZXZXZXZXZXZXY";
        pat = "XZXZXY";
        assertEquals(10, obj.KMPSearch(pat, txt), 0.0);
    }

    @Test
    public void testBinarySearchFound() {
        int[] arr = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        int target = 5;
        assertEquals(4, obj.binarySearch(arr, target));  // Target 5 is at index 4
    }

    @Test
    public void testBinarySearchNotFound() {
        int[] arr = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        int target = 11;
        assertEquals(-1, obj.binarySearch(arr, target));  // Target 11 is not in the array
    }

    @Test
    public void LinearTargetFound() {
        int[] array = {1, 2, 3, 4, 5};
        int target = 3;
        assertEquals(2, obj.linearSearch(array, target));
    }

    @Test
    public void LinearTargetNotFound() {
        int[] array = {10, 20, 30, 40};
        int target = 25;
        assertEquals(-1, obj.linearSearch(array, target));
    }

    @Test
    public void TestLCS() {
        String s1 = "abcdefghij";
        String s2 = "klmnopqrst";

        char[] X = s1.toCharArray();
        char[] Y = s2.toCharArray();
        int m = X.length;
        int n = Y.length;
        assertEquals(0, obj.LCS(X, Y, m, n), 0.0);

        s1 = "abcdefghijklmnopqrstuvwxyz";
        s2 = "zyxwvutsrqponmlkjihgfedcba";
        X = s1.toCharArray();
        Y = s2.toCharArray();
        m = X.length;
        n = Y.length;
        assertEquals(1, obj.LCS(X, Y, m, n), 0.0);

        s1 = "abcde";
        s2 = "abcdeabcde";
        X = s1.toCharArray();
        Y = s2.toCharArray();
        m = X.length;
        n = Y.length;
        assertEquals(5, obj.LCS(X, Y, m, n), 0.0);
    }

    @Test
    public void TestLongestPalindromicSubsequence() {
        String seq = "abcdefghijklmnopqrstuvwxy";
        int n = seq.length();
        assertEquals(1, obj.lps(seq.toCharArray(), 0, n - 1), 0.0);

        seq = "abcdefghijabcdefghijabcdefghij";
        n = seq.length();
        assertEquals(5, obj.lps(seq.toCharArray(), 0, n - 1), 0.0);

        seq = "a";
        n = seq.length();
        assertEquals(1, obj.lps(seq.toCharArray(), 0, n - 1), 0.0);

        seq = "abcdeffedcbazyxwvutsrqponmlkjihgfedcba";
        n = seq.length();
        assertEquals(13, obj.lps(seq.toCharArray(), 0, n - 1), 0.0);
    }

    @Test
    public void TestShortestCommonSequence() {
        String X = "AGGTB";
        String Y = "GXTXAYB";

        int[][] lookup = new int[X.length() + 1][Y.length() + 1];

        assertEquals(9, obj.superSeq(X, Y, X.length(), Y.length(), lookup), 0.0);

        X = "apqrstu";
        Y = "kplrmntuo";
        lookup = new int[X.length() + 1][Y.length() + 1];
        assertEquals(12, obj.superSeq(X, Y, X.length(), Y.length(), lookup), 0.0);
    }

    @Test
    public void TestmaxCommStr() {

        String s1 = "abcdef";
        String s2 = "zabcy";
        assertEquals(3, obj.maxCommStr(s1, s2), 0.0);

        String s3 = "abc";
        String s4 = "def";
        assertEquals(0, obj.maxCommStr(s3, s4), 0.0);
    }

    @Test
    public void TestLongestCommonPrefix() {
        String[] input = { "abcdefghijklmnopqrsxyz", "abcdefghijklmnopqrstuvxyz", "abcdefghijklmnopqrstuvwxyz" };
        assertEquals("abcdefghijklmnopqrs", obj.longestCommonPrefix(input));

        String[] input1 = { "", "abcdefghijklmno", "xyzabcdefghijklmnopqrs" };
        assertEquals("", obj.longestCommonPrefix(input1));

        String[] input2 = { "abcd", "efg", "zzzz" };
        assertEquals("", obj.longestCommonPrefix(input2));

        String[] input3 = { "abcdef" };
        assertEquals("abcdef", obj.longestCommonPrefix(input3));

        String[] input4 = { "helloworld", "helloworld!", "helloworldagain" };
        assertEquals("helloworld", obj.longestCommonPrefix(input4));

        String[] input5 = { "applepieistasty", "apricotsaregood", "apexpredator" };
        assertEquals("ap", obj.longestCommonPrefix(input5));
    }

    @Test
    public void TestLVP() {
        assertEquals(8, obj.findLongestValidParentheses("((()()()()(((())", 16), 0.0);
        assertEquals(0, obj.findLongestValidParentheses(")(", 2));
        assertEquals(0, obj.findLongestValidParentheses("", 0));
        assertEquals(0, obj.findLongestValidParentheses("))))))", 6), 0.0);
        assertEquals(4, obj.findLongestValidParentheses("((((())", 7), 0.0);
    }

    @Test
    public void TestEditDistance() {
        String str1 = "";
        String str2 = "";

        int n = str1.length(), m = str2.length();
        int[][] dp = new int[n + 1][m + 1];
        for (int i = 0; i < n + 1; i++)
            Arrays.fill(dp[i], -1);
        assertEquals(0, obj.calculateEditDistance(str1, str2, n, m, dp), 0.0);

        str1 = "algorithm";
        str2 = "altruistic";

        n = str1.length();
        m = str2.length();
        dp = new int[n + 1][m + 1];
        for (int i = 0; i < n + 1; i++)
            Arrays.fill(dp[i], -1);
        assertEquals(6, obj.calculateEditDistance(str1, str2, n, m, dp), 0.0);

        str1 = "editdistanceusingdynamicprogramming";
        str2 = "edit";

        n = str1.length();
        m = str2.length();
        dp = new int[n + 1][m + 1];
        for (int i = 0; i < n + 1; i++)
            Arrays.fill(dp[i], -1);
        assertEquals(31, obj.calculateEditDistance(str1, str2, n, m, dp), 0.0);

    }

    @Test
    public void TestLongestPalindromicSubstring() {
        String text = "ababaabacabacaba";
        assertEquals("abacabacaba", obj.findLongestPalindromicSubstring(text));
        text = "abcdeghghijkalmnahgyuwdggggaygufgsiugfisuuuuaaagggauuuuuuuhhhabbbbiiiiiiiiisssssssssoooooooaaauuuuuu";
        assertEquals("iiiiiiiii", obj.findLongestPalindromicSubstring(text));

        text = "abcddcbax";
        assertEquals("abcddcba", obj.findLongestPalindromicSubstring(text));

        text = "xyzabcdefghijklmnopqrsmadam";
        assertEquals("madam", obj.findLongestPalindromicSubstring(text));

    }

    @Test
    public void TestBoyerMoore() {
        // Test 1: Empty Text and Pattern
        char txt[] = "".toCharArray();
        char pat[] = "".toCharArray();
        assertEquals(0, obj.BoyerMoore(txt, pat));

        // Test 2: Pattern at the End
        txt = "ABCDEF".toCharArray();
        pat = "DEF".toCharArray();
        assertEquals(3, obj.BoyerMoore(txt, pat));

        // Test 3: Pattern Not Present
        txt = "ABCDEF".toCharArray();
        pat = "XYZ".toCharArray();
        assertEquals(-1, obj.BoyerMoore(txt, pat));

        // Test 4: Repeated Pattern
        txt = "ABABAB".toCharArray();
        pat = "AB".toCharArray();
        assertEquals(0, obj.BoyerMoore(txt, pat));
    }

    @Test
    public void TestSequenceAlignment() {
        String gene1 = "AGGGCT";
        String gene2 = "AGGCA";
        int misMatchPenalty = 3;
        int gapPenalty = 2;
        assertEquals(5, obj.SequenceAlignment(gene1, gene2, misMatchPenalty, gapPenalty));

        gene1 = "CG";
        gene2 = "CA";
        misMatchPenalty = 3;
        gapPenalty = 7;
        assertEquals(3, obj.SequenceAlignment(gene1, gene2, misMatchPenalty, gapPenalty));

        gene1 = "CG";
        gene2 = "CA";
        misMatchPenalty = 3;
        gapPenalty = 5;
        assertEquals(3, obj.SequenceAlignment(gene1, gene2, misMatchPenalty, gapPenalty));
    }

    @Test
    public void TestWildcardPattern() {
        String str = "baaabab";
        String pattern1 = "*****ba*****ab";
        assertEquals(true, obj.WildcardPattern(str, pattern1, str.length(), pattern1.length()));
        String pattern2 = "baaa?ab";
        assertEquals(true, obj.WildcardPattern(str, pattern2, str.length(), pattern2.length()));
        String pattern3 = "ba*a?";
        assertEquals(true, obj.WildcardPattern(str, pattern3, str.length(), pattern3.length()));
        String pattern4 = "a*ab";
        assertEquals(false, obj.WildcardPattern(str, pattern4, str.length(), pattern4.length()));
    }

    @Test
    public void TestminPalPartition() {
        String str = "geek";
        assertEquals(2, obj.minPalPartition(str));
        str = "aaaa";
        assertEquals(0, obj.minPalPartition(str));
        str = "abcde";
        assertEquals(4, obj.minPalPartition(str));
        str = "abbac";
        assertEquals(1, obj.minPalPartition(str));
    }

    @Test
    public void TestLongestRepeatingSubSeq() {
        String str = "aabb";
        assertEquals(2, obj.LongestRepeatingSubSeq(str));
        str = "abc";
        assertEquals(0, obj.LongestRepeatingSubSeq(str));
        str = "aab";
        assertEquals(1, obj.LongestRepeatingSubSeq(str));
        str = "axxxy";
        assertEquals(2, obj.LongestRepeatingSubSeq(str));
    }

    @Test
    public void TestLongestPrefixSuffix() {
        String str = "level";
        assertEquals(1, obj.longestPrefixSuffix(str));
        str = "aabcdaabc";
        assertEquals(4, obj.longestPrefixSuffix(str));
        str = "sgfgsdgdfsdd";
        assertEquals(0, obj.longestPrefixSuffix(str));
    }

    @Test
    public void TestKVowelWords() {
        int N = 1;
        int K = 0;
        assertEquals(21, obj.KVowelWords(N, K));
        N = 1;
        K = 1;
        assertEquals(26, obj.KVowelWords(N, K));
    }

    @Test
    public void TestLRR() {
        String str = "ncbjknsdjkcnsjkancjksdncjksdncjk";
        assertEquals("jknsdjkcnsjkancjksdncjksdncjkncb", obj.leftrotate(str, 3));
        assertEquals("jkncbjknsdjkcnsjkancjksdncjksdnc", obj.rightrotate(str, 2));
        str = "GeeksforGeeks";
        assertEquals("eksforGeeksGe", obj.leftrotate(str, 2));
        assertEquals("ksGeeksforGee", obj.rightrotate(str, 2));
        str = "qwertyu";
        assertEquals("ertyuqw", obj.leftrotate(str, 2));
        assertEquals("yuqwert", obj.rightrotate(str, 2));
    }

    @Test
    public void TestReverseVowel() {
        String str = "leetcode";
        assertEquals("leotcede", obj.reverseVowel(str));
        str = "aeiou";
        assertEquals("uoiea", obj.reverseVowel(str));
        str = "a";
        assertEquals("a", obj.reverseVowel(str));
    }

    @Test
    public void TestRepeatedStringMatch() {
        String a = "abcd", b = "cdabcdab";
        assertEquals(3, obj.repeatedStringMatch(a, b));
        a = "abcde";
        b = "bcdeabcdeabcdeabcdeabcd";
        assertEquals(5, obj.repeatedStringMatch(a, b));
        a = "cdmnslkd";
        b = "iwomeiofnmwo";
        assertEquals(-1, obj.repeatedStringMatch(a, b));
    }

    @Test
    public void testFindAllConcatenatedWordsInADict() {
        // Test Case 1: Empty array
        String[] words1 = {};
        assertEquals(0, obj.findAllConcatenatedWordsInADict(words1).size());

        // Test Case 2: No concatenated words
        String[] words2 = { "apple", "banana", "orange" };
        assertEquals(0, obj.findAllConcatenatedWordsInADict(words2).size());

        // Test Case 3: Single concatenated word
        String[] words3 = { "cat", "dog", "catdog" };
        assertEquals(1, obj.findAllConcatenatedWordsInADict(words3).size());
        assertTrue(obj.findAllConcatenatedWordsInADict(words3).contains("catdog"));

        // Test Case 4: Multiple concatenated words
        String[] words4 = { "hello", "world", "helloworld", "concat", "enated", "concatenated" };
        assertEquals(2, obj.findAllConcatenatedWordsInADict(words4).size());
        assertTrue(obj.findAllConcatenatedWordsInADict(words4).contains("helloworld"));
        assertTrue(obj.findAllConcatenatedWordsInADict(words4).contains("concatenated"));

        // Test Case 5: Empty string in the array
        String[] words5 = { "cat", "dog", "catdog" };
        assertEquals(1, obj.findAllConcatenatedWordsInADict(words5).size());
        assertTrue(obj.findAllConcatenatedWordsInADict(words5).contains("catdog"));

        // Test Case 6: Words with spaces
        String[] words6 = { "cat", "dog", "cat dog", "dog cat", "catdog" };
        assertEquals(1, obj.findAllConcatenatedWordsInADict(words6).size());
        assertTrue(obj.findAllConcatenatedWordsInADict(words6).contains("catdog"));

        // Test Case 7: Words with special characters
        String[] words7 = { "cat", "dog", "cat@dog", "dog@cat", "catdog" };
        assertEquals(1, obj.findAllConcatenatedWordsInADict(words7).size());
        assertTrue(obj.findAllConcatenatedWordsInADict(words7).contains("catdog"));
    }

    @Test
    public void testWordBreak1() {
        // Test Case 2: No valid word break
        String s2 = "catsandog";
        List<String> wordDict2 = Arrays.asList("cats", "dog", "sand", "and", "cat");
        assertFalse(obj.wordBreak1(s2, wordDict2));

        // Test Case 3: Valid word break
        String s3 = "applepenapple";
        List<String> wordDict3 = Arrays.asList("apple", "pen");
        assertTrue(obj.wordBreak1(s3, wordDict3));

        // Test Case 4: Complex word break
        String s4 = "catsanddog";
        List<String> wordDict4 = Arrays.asList("cats", "dog", "sand", "and", "cat");
        assertTrue(obj.wordBreak1(s4, wordDict4));

        // Test Case 5: Word break with trailing and leading spaces
        String s5 = "applepenapple";
        List<String> wordDict5 = Arrays.asList("apple", "pen");
        assertTrue(obj.wordBreak1(s5, wordDict5));

        // Test Case 6: Empty string and empty dictionary
        String s6 = "";
        List<String> wordDict6 = Collections.emptyList();
        assertTrue(obj.wordBreak1(s6, wordDict6));
    }

    @Test
    public void testWordBreak2() {
        // Test Case 2: No valid word break
        String s2 = "catsandog";
        List<String> wordDict2 = Arrays.asList("cats", "dog", "sand", "and", "cat");
        assertEquals(0, obj.wordBreak2(s2, wordDict2).size());

        // Test Case 3: Valid word break
        String s3 = "catsanddog";
        List<String> wordDict3 = Arrays.asList("cat", "cats", "and", "sand", "dog");
        assertEquals(Arrays.asList("cat sand dog", "cats and dog"), obj.wordBreak2(s3, wordDict3));

        // Test Case 4: Complex word break with duplicates
        String s4 = "aaa";
        List<String> wordDict4 = Arrays.asList("a", "aa");
        assertEquals(Arrays.asList("a aa", "a a a", "aa a"), obj.wordBreak2(s4, wordDict4));

        // Test Case 5: Word break with trailing and leading spaces
        String s5 = "catsanddog";
        List<String> wordDict5 = Arrays.asList("cat", "cats", "and", "sand", "dog");
        assertEquals(Arrays.asList("cat sand dog", "cats and dog"), obj.wordBreak2(s5, wordDict5));
    }

    @Test
    public void testAtMostNGivenDigitSet() {
        // Test Case 1
        String[] digits1 = { "1", "3", "5", "7" };
        int n1 = 100;
        assertEquals(20, obj.atMostNGivenDigitSet(digits1, n1));

        // Test Case 2
        String[] digits2 = { "1", "4", "9" };
        int n2 = 1000000000;
        assertEquals(29523, obj.atMostNGivenDigitSet(digits2, n2));

        // Test Case 3
        String[] digits3 = { "7" };
        int n3 = 8;
        assertEquals(1, obj.atMostNGivenDigitSet(digits3, n3));

    }

    @Test
    public void TestpalindromePairs() {
        String words[] = { "abcd", "dcba", "lls", "s", "sssll" };
        List<List<Integer>> ans = new ArrayList<>(Arrays.asList(
                Arrays.asList(0, 1),
                Arrays.asList(1, 0),
                Arrays.asList(3, 2),
                Arrays.asList(2, 4)));
        assertEquals(ans, obj.palindromePairs(words));
        String words1[] = { "bat", "tab", "cat" };
        List<List<Integer>> ans1 = new ArrayList<>(Arrays.asList(
                Arrays.asList(0, 1),
                Arrays.asList(1, 0)));
        assertEquals(ans1, obj.palindromePairs(words1));
        String words2[] = { "a", "" };
        List<List<Integer>> ans2 = new ArrayList<>(Arrays.asList(
                Arrays.asList(1, 0),
                Arrays.asList(0, 1)));
        assertEquals(ans2, obj.palindromePairs(words2));
    }

    @Test
    public void TestminStickers() {
        String stickers[] = { "with", "example", "science" };
        String target = "thehat";
        assertEquals(3, obj.minStickers(stickers, target));
        String stickers1[] = { "notice", "possible" };
        String target1 = "basicbasic";
        assertEquals(-1, obj.minStickers(stickers1, target1));
    }

    @Test
    public void TestkSimilarity() {
        String a = "ab", b = "ba";
        assertEquals(1, obj.kSimilarity(a, b));
        a = "abc";
        b = "bca";
        assertEquals(2, obj.kSimilarity(a, b));
    }
}
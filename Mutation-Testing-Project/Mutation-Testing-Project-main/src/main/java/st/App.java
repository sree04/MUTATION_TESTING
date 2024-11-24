package st;

import java.util.*;

public class App {

	int rabinKarp(String pattern, String text, int prime) {
		int base = 256; // Number of possible characters in the input
		int patternLength = pattern.length();
		int textLength = text.length();
		int patternHash = 0; // Hash value for the pattern
		int windowHash = 0;
		int hashFactor = 1;
		for (int i = 0; i < patternLength - 1; i++) {
			hashFactor = (hashFactor * base) % prime;
		}

		for (int i = 0; i < patternLength; i++) {
			patternHash = (base * patternHash + pattern.charAt(i)) % prime;
			windowHash = (base * windowHash + text.charAt(i)) % prime;
		}

		for (int i = 0; i <= textLength - patternLength; i++) {
			if (patternHash == windowHash) {
				boolean match = true;
				for (int j = 0; j < patternLength; j++) {
					if (text.charAt(i + j) != pattern.charAt(j)) {
						match = false;
						break;
					}
				}
				if (match) {
					return i;
				}
			}

			if (i < textLength - patternLength) {
				windowHash = (base * (windowHash - text.charAt(i) * hashFactor) + text.charAt(i + patternLength)) % prime;


				if (windowHash < 0) {
					windowHash += prime;
				}
			}
		}
		return -1;
	}


	int ZAlgorithm(String text, String pattern) {
		// Create concatenated string "P$T"
		String concat = pattern + "$" + text;
		int l = concat.length();
		int[] Z = new int[l];

		// Initialize the Z array
		int L = 0, R = 0;

		// Calculate Z array for the concatenated string
		for (int i = 1; i < l; i++) {
			if (i > R) {
				// If we're outside the current [L, R] range, start fresh
				L = R = i;
				while (R < l && concat.charAt(R) == concat.charAt(R - L)) {
					R++;
				}
				Z[i] = R - L;
				R--;
			} else {
				int k = i - L;
				if (Z[k] < R - i + 1) {
					// Z[k] is valid within the range [L, R]
					Z[i] = Z[k];
				} else {
					// Need to manually extend the match beyond R
					L = i;
					while (R < l && concat.charAt(R) == concat.charAt(R - L)) {
						R++;
					}
					Z[i] = R - L;
					R--;
				}
			}
		}

		// Now search for the pattern's occurrence using the Z array
		for (int i = 0; i < l; i++) {
			if (Z[i] == pattern.length()) {
				// Match found at position i - pattern.length() - 1 in the original text
				return i - pattern.length() - 1;
			}
		}

		// If no match is found, return -1
		return -1;
	}

	int binarySearch(int[] arr, int target) {
		int left = 0;
		int right = arr.length - 1;

		while (left <= right) {
			int mid = left + (right - left) / 2;
			if (arr[mid] == target) {
				return mid;
			}
			if (arr[mid] > target) {
				right = mid - 1;
			}
			else {
				left = mid + 1;
			}
		}

		// Target not found
		return -1;
	}

	int linearSearch(int[] arr, int target) {
		for (int i = 0; i < arr.length; i++) {
			if (arr[i] == target) {
				return i;
			}
		}
		return -1;
	}

	public int KMPSearch(String pattern, String text) {
		int patternLength = pattern.length();
		int textLength = text.length();

		// Compute the LPS (Longest Prefix Suffix) array
		int[] lps = new int[patternLength];
		int prefixLength = 0;
		lps[0] = 0;

		// Build the LPS array
		for (int i = 1; i < patternLength; ) {
			if (pattern.charAt(i) == pattern.charAt(prefixLength)) {
				prefixLength++;
				lps[i] = prefixLength;
				i++;
			} else {
				if (prefixLength != 0) {
					prefixLength = lps[prefixLength - 1];
				} else {
					lps[i] = 0;
					i++;
				}
			}
		}

		// Search for the pattern in the text
		int textIndex = 0, patternIndex = 0;
		while (textIndex < textLength) {
			if (pattern.charAt(patternIndex) == text.charAt(textIndex)) {
				textIndex++;
				patternIndex++;
			}

			if (patternIndex == patternLength) {
				return textIndex - patternIndex; // Match found
			} else if (textIndex < textLength && pattern.charAt(patternIndex) != text.charAt(textIndex)) {
				if (patternIndex != 0) {
					patternIndex = lps[patternIndex - 1];
				} else {
					textIndex++;
				}
			}
		}

		return -1; // No match found
	}


	int LCS(char[] X, char[] Y, int m, int n) {
		int L[][] = new int[m + 1][n + 1];

		/*
		 * Following steps build L[m+1][n+1] in bottom up fashion. Note
		 * that L[i][j] contains length of LCS of X[0..i-1] and Y[0..j-1]
		 */
		for (int i = 0; i <= m; i++) {
			for (int j = 0; j <= n; j++) {
				if (i == 0 || j == 0)
					L[i][j] = 0;
				else if (X[i - 1] == Y[j - 1])
					L[i][j] = L[i - 1][j - 1] + 1;
				else
					L[i][j] = max(L[i - 1][j], L[i][j - 1]);
			}
		}
		return L[m][n];
	}

	/* Utility function to get max of 2 integers */
	int max(int a, int b) {
		return (a > b) ? a : b;
	}

	int lps(char seq[], int i, int j) {
		// Base Case 1: If there is only 1 character
		if (i == j) {
			return 1;
		}

		// Base Case 2: If there are only 2 characters and both are same
		if (seq[i] == seq[j] && i + 1 == j) {
			return 2;
		}

		// If the first and last characters match
		if (seq[i] == seq[j]) {
			return lps(seq, i + 1, j - 1) + 2;
		}

		// If the first and last characters do not match
		return max(lps(seq, i, j - 1), lps(seq, i + 1, j));
	}

	int superSeq(String X, String Y, int n, int m, int[][] lookup) {

		if (m == 0 || n == 0) {
			lookup[n][m] = n + m;
		}

		if (lookup[n][m] == 0)
			if (X.charAt(n - 1) == Y.charAt(m - 1)) {
				lookup[n][m] = superSeq(X, Y, n - 1, m - 1, lookup)
						+ 1;
			}

			else {
				lookup[n][m] = Math.min(
						superSeq(X, Y, n - 1, m, lookup) + 1,
						superSeq(X, Y, n, m - 1, lookup) + 1);
			}

		return lookup[n][m];
	}

	int maxCommStr(String s1, String s2) {
		int m = s1.length();
		int n = s2.length();
		int[][] dp = new int[m + 1][n + 1];
		int maxLength = 0;

		for (int i = 1; i <= m; i++) {
			for (int j = 1; j <= n; j++) {
				if (s1.charAt(i - 1) == s2.charAt(j - 1)) {
					dp[i][j] = dp[i - 1][j - 1] + 1;
					maxLength = Math.max(maxLength, dp[i][j]);
				}
			}
		}

		return maxLength;
	}
	String longestCommonPrefix(String[] strings) {
		int length = strings.length;

		// Return an empty string if the array is empty
		if (length == 0) {
			return "";
		}

		// If the array has only one string, return it
		if (length == 1) {
			return strings[0];
		}

		// Sort the array of strings to bring lexicographically similar strings closer
		Arrays.sort(strings);

		// Find the smallest length between the first and last strings in the sorted array
		int minLength = Math.min(strings[0].length(), strings[length - 1].length());

		// Compare characters of the first and last strings up to the minimum length
		int index = 0;
		while (index < minLength && strings[0].charAt(index) == strings[length - 1].charAt(index)) {
			index++;
		}

		// Extract the substring from the first string up to the matched index
		return strings[0].substring(0, index);
	}

	int findLongestValidParentheses(String str, int length) {
		// Track counts of opening and closing parentheses
		int openCount = 0, closeCount = 0;
		int maxLength = 0;

		// Scan the string from left to right
		for (int index = 0; index < length; index++) {
			// Update counters based on the character
			if (str.charAt(index) == '(') {
				openCount++;
			} else {
				closeCount++;
			}

			// When counts are equal, update the maximum length
			if (openCount == closeCount) {
				maxLength = Math.max(maxLength, 2 * closeCount);
			}
			// Reset counters if the sequence becomes invalid
			else if (closeCount > openCount) {
				openCount = closeCount = 0;
			}
		}

		// Reset counters for the reverse traversal
		openCount = closeCount = 0;

		// Scan the string from right to left
		for (int index = length - 1; index >= 0; index--) {
			// Update counters based on the character
			if (str.charAt(index) == '(') {
				openCount++;
			} else {
				closeCount++;
			}

			// When counts are equal, update the maximum length
			if (openCount == closeCount) {
				maxLength = Math.max(maxLength, 2 * openCount);
			}
			// Reset counters if the sequence becomes invalid
			else if (openCount > closeCount) {
				openCount = closeCount = 0;
			}
		}

		return maxLength;
	}

	public int calculateEditDistance(String str1, String str2, int len1, int len2, int[][] memo) {
		// Base cases: if one string is empty, return the length of the other
		if (len1 == 0) {
			return len2;
		}
		if (len2 == 0) {
			return len1;
		}

		// Check if the result for this state has already been computed
		if (memo[len1][len2] != -1) {
			return memo[len1][len2];
		}

		// If the characters are the same, move to the next pair
		if (str1.charAt(len1 - 1) == str2.charAt(len2 - 1)) {
			memo[len1][len2] = calculateEditDistance(str1, str2, len1 - 1, len2 - 1, memo);
		} else {
			// Calculate the cost for each possible operation
			int costInsert = calculateEditDistance(str1, str2, len1, len2 - 1, memo);
			int costDelete = calculateEditDistance(str1, str2, len1 - 1, len2, memo);
			int costReplace = calculateEditDistance(str1, str2, len1 - 1, len2 - 1, memo);

			// Choose the operation with the minimum cost and add 1 for the current step
			memo[len1][len2] = 1 + Math.min(costInsert, Math.min(costDelete, costReplace));
		}

		return memo[len1][len2];
	}

	public String findLongestPalindromicSubstring(String input) {
		int length = input.length();
		if (length == 0) {
			return null;
		}

		// Transform the string to accommodate palindromes of both even and odd lengths
		int transformedLength = 2 * length + 1;
		int[] lpsArray = new int[transformedLength];
		lpsArray[0] = 0;
		lpsArray[1] = 1;

		int currentCenter = 1; // Center of the current palindrome
		int rightBoundary = 2; // Right boundary of the current palindrome
		int maxLPSLength = 0;
		int maxLPSCenter = 0;

		for (int pos = 2; pos < transformedLength; pos++) {
			// Calculate the mirrored position of pos with respect to the current center
			int mirrorPos = 2 * currentCenter - pos;
			int boundaryDifference = rightBoundary - pos;

			// Use previously computed palindrome lengths to initialize the current length
			if (boundaryDifference > 0) {
				lpsArray[pos] = Math.min(lpsArray[mirrorPos], boundaryDifference);
			} else {
				lpsArray[pos] = 0;
			}

			// Attempt to expand the palindrome centered at pos
			while ((pos + lpsArray[pos] + 1 < transformedLength && pos - lpsArray[pos] - 1 >= 0) &&
					((pos + lpsArray[pos] + 1) % 2 == 0 ||
							input.charAt((pos + lpsArray[pos] + 1) / 2) == input.charAt((pos - lpsArray[pos] - 1) / 2))) {
				lpsArray[pos]++;
			}

			// Update the maximum LPS length and center if needed
			if (lpsArray[pos] > maxLPSLength) {
				maxLPSLength = lpsArray[pos];
				maxLPSCenter = pos;
			}

			// Adjust the center and right boundary if the palindrome at pos extends beyond the current right boundary
			if (pos + lpsArray[pos] > rightBoundary) {
				currentCenter = pos;
				rightBoundary = pos + lpsArray[pos];
			}
		}

		// Calculate the start and end indices of the longest palindromic substring
		int start = (maxLPSCenter - maxLPSLength) / 2;
		int end = start + maxLPSLength - 1;

		return input.substring(start, end + 1);
	}


	// BOYER MOORE ALGORITHM
	int NO_OF_CHARS = 256;

	// The preprocessing function for Boyer Moore's
	// bad character heuristic
	void badCharHeuristic(char[] str, int size, int badchar[]) {

		// Initialize all occurrences as -1
		for (int i = 0; i < NO_OF_CHARS; i++)
			badchar[i] = -1;

		// Fill the actual value of last occurrence
		// of a character (indices of table are ascii and values are index of
		// occurrence)
		for (int i = 0; i < size; i++)
			badchar[(int) str[i]] = i;
	}

	/*
	 * A pattern searching function that uses Bad
	 * Character Heuristic of Boyer Moore Algorithm
	 */
	int BoyerMoore(char txt[], char pat[]) {
		int m = pat.length;
		int n = txt.length;

		int badchar[] = new int[NO_OF_CHARS];

		/*
		 * Fill the bad character array by calling
		 * the preprocessing function badCharHeuristic()
		 * for given pattern
		 */
		badCharHeuristic(pat, m, badchar);

		int s = 0; // s is shift of the pattern with
					// respect to text
		// there are n-m+1 potential alignments
		while (s <= (n - m)) {
			int j = m - 1;

			/*
			 * Keep reducing index j of pattern while
			 * characters of pattern and text are
			 * matching at this shift s
			 */
			while (j >= 0 && pat[j] == txt[s + j])
				j--;

			/*
			 * If the pattern is present at current
			 * shift, then index j will become -1 after
			 * the above loop
			 */
			if (j < 0) {
				return s;

				/*
				 * Shift the pattern so that the next
				 * character in text aligns with the last
				 * occurrence of it in pattern.
				 * The condition s+m < n is necessary for
				 * the case when pattern occurs at the end
				 * of text
				 */
				// txt[s+m] is character after the pattern in text
				// s += (s+m < n)? m-badchar[txt[s+m]] : 1;

			} else
				/*
				 * Shift the pattern so that the bad character
				 * in text aligns with the last occurrence of
				 * it in pattern. The max function is used to
				 * make sure that we get a positive shift.
				 * We may get a negative shift if the last
				 * occurrence of bad character in pattern
				 * is on the right side of the current
				 * character.
				 */
				s += max(1, j - badchar[txt[s + j]]);
		}
		return -1;
	}

	int SequenceAlignment(String x, String y, int pxy, int pgap) {
		int i, j; // initialising variables

		int m = x.length(); // length of gene1
		int n = y.length(); // length of gene2

		// table for storing optimal
		// substructure answers
		int dp[][] = new int[n + m + 1][n + m + 1];

		for (int[] x1 : dp)
			Arrays.fill(x1, 0);

		// initialising the table
		for (i = 0; i <= (n + m); i++) {
			dp[i][0] = i * pgap;
			dp[0][i] = i * pgap;
		}

		// calculating the
		// minimum penalty
		for (i = 1; i <= m; i++) {
			for (j = 1; j <= n; j++) {
				if (x.charAt(i - 1) == y.charAt(j - 1)) {
					dp[i][j] = dp[i - 1][j - 1];
				} else {
					dp[i][j] = Math.min(Math.min(dp[i - 1][j - 1] + pxy,
							dp[i - 1][j] + pgap),
							dp[i][j - 1] + pgap);
				}
			}
		}

		// Reconstructing the solution
		int l = n + m; // maximum possible length

		i = m;
		j = n;

		int xpos = l;
		int ypos = l;

		// Final answers for
		// the respective strings
		int xans[] = new int[l + 1];
		int yans[] = new int[l + 1];

		while (!(i == 0 || j == 0)) {
			if (x.charAt(i - 1) == y.charAt(j - 1)) {
				xans[xpos--] = (int) x.charAt(i - 1);
				yans[ypos--] = (int) y.charAt(j - 1);
				i--;
				j--;
			} else if (dp[i - 1][j - 1] + pxy == dp[i][j]) {
				xans[xpos--] = (int) x.charAt(i - 1);
				yans[ypos--] = (int) y.charAt(j - 1);
				i--;
				j--;
			} else if (dp[i - 1][j] + pgap == dp[i][j]) {
				xans[xpos--] = (int) x.charAt(i - 1);
				yans[ypos--] = (int) '_';
				i--;
			} else if (dp[i][j - 1] + pgap == dp[i][j]) {
				xans[xpos--] = (int) '_';
				yans[ypos--] = (int) y.charAt(j - 1);
				j--;
			}
		}
		while (xpos > 0) {
			if (i > 0)
				xans[xpos--] = (int) x.charAt(--i);
			else
				xans[xpos--] = (int) '_';
		}
		while (ypos > 0) {
			if (j > 0)
				yans[ypos--] = (int) y.charAt(--j);
			else
				yans[ypos--] = (int) '_';
		}

		// Since we have assumed the
		// answer to be n+m long,
		// we need to remove the extra
		// gaps in the starting id
		// represents the index from
		// which the arrays xans,
		// yans are useful
		int id = 1;
		for (i = l; i >= 1; i--) {
			if ((char) yans[i] == '_' && (char) xans[i] == '_') {
				id = i + 1;
				break;
			}
		}
		return dp[m][n];
	}

	boolean WildcardPattern(String str, String pattern,
			int n, int m) {
		// empty pattern can only match with
		// empty string
		if (m == 0)
			return (n == 0);

		// lookup table for storing results of
		// subproblems
		boolean[][] lookup = new boolean[n + 1][m + 1];

		// initialize lookup table to false
		for (int i = 0; i < n + 1; i++)
			Arrays.fill(lookup[i], false);

		// empty pattern can match with empty string
		lookup[0][0] = true;

		// Only '*' can match with empty string
		for (int j = 1; j <= m; j++)
			if (pattern.charAt(j - 1) == '*')
				lookup[0][j] = lookup[0][j - 1];

		// fill the table in bottom-up fashion
		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= m; j++) {
				// Two cases if we see a '*'
				// a) We ignore '*'' character and move
				// to next character in the pattern,
				// i.e., '*' indicates an empty
				// sequence.
				// b) '*' character matches with ith
				// character in input
				if (pattern.charAt(j - 1) == '*')
					lookup[i][j] = lookup[i][j - 1]
							|| lookup[i - 1][j];

				// Current characters are considered as
				// matching in two cases
				// (a) current character of pattern is '?'
				// (b) characters actually match
				else if (pattern.charAt(j - 1) == '?'
						|| str.charAt(i - 1) == pattern.charAt(j - 1))
					lookup[i][j] = lookup[i - 1][j - 1];

				// If characters don't match
				else
					lookup[i][j] = false;
			}
		}

		return lookup[n][m];
	}

	int minPalPartition(String str) {
		// Get the length of the string
		int n = str.length();

		/*
		 * Create two arrays to build the solution
		 * in bottom up manner
		 * C[i] = Minimum number of cuts needed for
		 * palindrome partitioning of substring
		 * str[0..i]
		 * P[i][j] = true if substring str[i..j] is
		 * palindrome, else false
		 * Note that C[i] is 0 if P[0][i] is true
		 */
		int[] C = new int[n];
		boolean[][] P = new boolean[n][n];

		int i, j, k, L; // different looping variables

		// Every substring of length 1 is a palindrome
		for (i = 0; i < n; i++) {
			P[i][i] = true;
		}

		/*
		 * L is substring length. Build the solution
		 * in bottom up manner by considering all substrings
		 * of length starting from 2 to n.
		 */
		for (L = 2; L <= n; L++) {
			// For substring of length L, set different
			// possible starting indexes
			for (i = 0; i < n - L + 1; i++) {
				j = i + L - 1; // Set ending index

				// If L is 2, then we just need to
				// compare two characters. Else need to
				// check two corner characters and value
				// of P[i+1][j-1]
				if (L == 2)
					P[i][j] = (str.charAt(i) == str.charAt(j));
				else
					P[i][j] = (str.charAt(i) == str.charAt(j)) && P[i + 1][j - 1];
			}
		}

		for (i = 0; i < n; i++) {
			if (P[0][i] == true)
				C[i] = 0;
			else {
				C[i] = Integer.MAX_VALUE;
				for (j = 0; j < i; j++) {
					if (P[j + 1][i] == true && 1 + C[j] < C[i])
						C[i] = 1 + C[j];
				}
			}
		}

		// Return the min cut value for complete
		// string. i.e., str[0..n-1]
		return C[n - 1];
	}

	int LongestRepeatingSubSeq(String str) {
		int n = str.length();

		// Create and initialize DP table
		int[][] dp = new int[n + 1][n + 1];

		// Fill dp table (similar to LCS loops)
		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= n; j++) {
				// If characters match and indexes are not same
				if (str.charAt(i - 1) == str.charAt(j - 1) && i != j)
					dp[i][j] = 1 + dp[i - 1][j - 1];

				// If characters do not match
				else
					dp[i][j] = Math.max(dp[i][j - 1], dp[i - 1][j]);
			}
		}
		return dp[n][n];
	}

	int longestPrefixSuffix(String s) {
		int n = s.length();

		int lps[] = new int[n];

		// lps[0] is always 0
		lps[0] = 0;

		// length of the previous
		// longest prefix suffix
		int len = 0;

		// the loop calculates lps[i]
		// for i = 1 to n-1
		int i = 1;
		while (i < n) {
			if (s.charAt(i) == s.charAt(len)) {
				len++;
				lps[i] = len;
				i++;
			}

			// (pat[i] != pat[len])
			else {
				// This is tricky. Consider
				// the example. AAACAAAA
				// and i = 7. The idea is
				// similar to search step.
				if (len != 0) {
					len = lps[len - 1];

					// Also, note that we do
					// not increment i here
				}

				// if (len == 0)
				else {
					lps[i] = 0;
					i++;
				}
			}
		}

		int res = lps[n - 1];

		// Since we are looking for
		// non overlapping parts.
		return (res > n / 2) ? n / 2 : res;
	}

	// Number of distinct words of size N with at most K contiguous vowels
	int power(int x, int y, int p) {
		int res = 1;
		x = x % p;

		if (x == 0)
			return 0;

		while (y > 0) {
			if ((y & 1) != 0)
				res = (res * x) % p;

			y = y >> 1;
			x = (x * x) % p;
		}
		return res;
	}

	int KVowelWords(int N, int K) {
		int i, j;
		int MOD = 1000000007;

		// Array dp to store number of ways
		int[][] dp = new int[N + 1][K + 1];

		int sum = 1;
		for (i = 1; i <= N; i++) {

			// dp[i][0] = (dp[i-1][0]+dp[i-1][1]..dp[i-1][k])*21
			dp[i][0] = sum * 21;
			dp[i][0] %= MOD;

			// Now setting sum to be dp[i][0]
			sum = dp[i][0];

			for (j = 1; j <= K; j++) {

				// If j>i, no ways are possible to create
				// a string with length i and vowel j
				if (j > i)
					dp[i][j] = 0;

				else if (j == i) {

					// If j = i all the character should
					// be vowel
					dp[i][j] = power(5, i, MOD);
				} else {

					// dp[i][j] relation with dp[i-1][j-1]
					dp[i][j] = dp[i - 1][j - 1] * 5;
				}

				dp[i][j] %= MOD;

				// Adding dp[i][j] in the sum
				sum += dp[i][j];
				sum %= MOD;
			}
		}
		return sum;
	}

	String leftrotate(String str1, int n) {

		// creating extended string and index for new
		// rotated string
		String temp = str1 + str1;
		int l1 = str1.length();

		String Lfirst = temp.substring(n, n + l1);

		// now returning string
		return Lfirst;
	}

	String rightrotate(String str1, int n) {
		return leftrotate(str1, str1.length() - n);
	}
	// LEFT AND RIGHT ROTATION OF A STRING

	boolean isVowel(char c) {
		return (c == 'a' || c == 'A' || c == 'e'
				|| c == 'E' || c == 'i' || c == 'I'
				|| c == 'o' || c == 'O' || c == 'u'
				|| c == 'U');
	}

	String reverseVowel(String str) {
		// Start two indexes from two corners
		// and move toward each other
		int i = 0;
		int j = str.length() - 1;
		char[] str1 = str.toCharArray();
		while (i < j) {
			if (!isVowel(str1[i])) {
				i++;
				continue;
			}
			if (!isVowel(str1[j])) {
				j--;
				continue;
			}

			// swapping
			char t = str1[i];
			str1[i] = str1[j];
			str1[j] = t;

			i++;
			j--;
		}
		String str2 = String.copyValueOf(str1);
		return str2;
	}

	// Horspool BM algorithm pattern searching
	int repeatedStringMatch(String a, String b) {
		StringBuilder text = new StringBuilder(a);
		char[] pat = b.toCharArray();
		int n = b.length();
		if (n == 0)
			return 0;
		int rep = 1;

		while (text.length() < n) {
			text.append(a);
			rep++;
		}

		int m = text.length();
		char[] ar = text.toString().toCharArray();

		int[] shifts = getShifts(b.toCharArray(), n);

		int i = n - 1;

		int contendor = horspool(ar, m, n, rep, shifts, pat);
		if (contendor > 0)
			return contendor;

		text.append(a);
		ar = (text.toString()).toCharArray();
		rep++;
		m = text.length();
		contendor = horspool(ar, m, n, rep, shifts, pat);
		if (contendor > 0)
			return contendor;

		return -1;
	}

	int horspool(char[] ar, int m, int n, int rep, int[] shifts, char[] pat) {
		int i = n - 1;

		while (i < m) {
			int k = i;
			int j = n - 1;

			i += shifts[ar[i] - 97];
			while (k >= 0 && j >= 0 && ar[k--] == pat[j])
				j--;

			if (j < 0)
				return rep;
		}
		return -1;
	}

	int[] getShifts(char[] ar, int m) {
		int[] res = new int[26];
		Arrays.fill(res, m);
		for (int i = 0; i < m - 1; i++)
			res[ar[i] - 97] = m - i - 1;

		return res;
	}

	private boolean dfs(final String word, int length, final boolean[] visited, final Set<String> dictionary) {
		if (length == word.length()) {
			return true;
		}

		if (visited[length]) {
			return false;
		}

		visited[length] = true;

		int end = (length == 0) ? word.length() - 1 : word.length();
		for (int i = end; i > length; --i) {
			String currentSubstring = word.substring(length, i);
			if (dictionary.contains(currentSubstring) && dfs(word, i, visited, dictionary)) {
				return true;
			}
		}

		return false;
	}

	public List<String> findAllConcatenatedWordsInADict(String[] words) {
		final Set<String> dictionary = new HashSet<>(Arrays.asList(words));
		final List<String> answer = new ArrayList<>();

		for (final String word : words) {
			final int length = word.length();
			final boolean[] visited = new boolean[length];

			if (dfs(word, 0, visited, dictionary)) {
				answer.add(word);
			}
		}

		return answer;
	}

	public boolean wordBreak1(String s, List<String> wordDict) {
		Map<String, Boolean> dp = new HashMap<>();
		Set<String> dict = new HashSet<>(wordDict);
		return util1(0, s.length() - 1, s, dict, dp);
	}

	public boolean util1(int i, int j, String s, Set<String> dict, Map<String, Boolean> dp) {
		if (i > j) {
			return true;
		}

		String key = i + "|" + j;
		if (dp.containsKey(key)) {
			return dp.get(key);
		}

		if (dict.contains(s.substring(i, j + 1))) {
			dp.put(key, true);
			return true;
		}

		boolean ret = false;
		for (int br = i; br <= j - 1; br++) {
			boolean cur = util1(i, br, s, dict, dp) && util1(br + 1, j, s, dict, dp);
			ret = ret || cur;
		}

		dp.put(key, ret);
		return ret;
	}

	public List<String> wordBreak2(String s, List<String> wordDict) {
		Map<String, List<String>> dp = new HashMap<>();
		Set<String> dict = new HashSet<>(wordDict);

		return util2(0, s.length() - 1, s, dict, dp);
	}

	public List<String> util2(int i, int j, String s, Set<String> dict, Map<String, List<String>> dp) {
		if (i > j) {
			List<String> ret = new ArrayList<>();
			ret.add(""); // Add an empty string to indicate a valid break
			return ret;
		}

		String key = i + "|" + j;
		if (dp.containsKey(key)) {
			return dp.get(key);
		}

		Set<String> retList = new HashSet<>();

		for (int br = i; br <= j - 1; br++) {
			List<String> left = util2(i, br, s, dict, dp);
			List<String> right = util2(br + 1, j, s, dict, dp);

			if (!left.isEmpty() && !right.isEmpty()) {
				for (String l : left) {
					for (String r : right) {
						String toAdd = l + " " + r;
						retList.add(toAdd);
					}
				}
			}
		}

		List<String> ret = new ArrayList<>(retList);

		// Check if the substring from i to j is a valid word in the dictionary
		if (dict.contains(s.substring(i, j + 1))) {
			ret.add(s.substring(i, j + 1));
		}

		dp.put(key, ret);
		return ret;
	}

	public int atMostNGivenDigitSet(String[] D, int N) {
		int k = D.length;
		int[] digits = new int[k];
		for (int i = 0; i < k; i++) {
			digits[i] = Integer.parseInt(D[i]);
		}

		int cnt = 0;
		String s = String.valueOf(N);
		int len = s.length();
		int[] rates = new int[len];
		rates[0] = 1;

		// Calculate rates for less digits: k + k^2 + .. + k^(len-1)
		for (int i = 0; i < len - 1; i++) {
			rates[i + 1] = rates[i] * k;
			cnt += rates[i + 1];
		}

		// Add count for same digits
		cnt += helper(digits, rates, s);
		return cnt;
	}

	private int helper(int[] digits, int[] rates, String s) {
		if (s.length() == 0) {
			return 1;
		}

		int n = s.charAt(0) - '0';
		int cnt = countLT(digits, n) * rates[s.length() - 1];

		if (exists(digits, n)) {
			cnt += helper(digits, rates, s.substring(1));
		}

		return cnt;
	}

	private int countLT(int[] digits, int n) {
		int cnt = 0;
		for (int d : digits) {
			if (d < n) {
				cnt++;
			} else {
				break;
			}
		}
		return cnt;
	}

	private boolean exists(int[] digits, int n) {
		for (int d : digits) {
			if (d == n) {
				return true;
			}
			if (d > n) {
				break;
			}
		}
		return false;
	}

	public String reverseStr(String str) {
		StringBuilder sb = new StringBuilder(str);
		return sb.reverse().toString();
	}

	public boolean isPalindrome(String s) {
		int i = 0;
		int j = s.length() - 1;
		while (i <= j) {
			if (s.charAt(i) != s.charAt(j)) {
				return false;
			}
			i++;
			j--;
		}
		return true;
	}

	public List<List<Integer>> palindromePairs(String[] words) {
		List<List<Integer>> res = new ArrayList<List<Integer>>();
		if (words == null || words.length == 0) {
			return res;
		}
		// build the map save the key-val pairs: String - idx
		HashMap<String, Integer> map = new HashMap<>();
		for (int i = 0; i < words.length; i++) {
			map.put(words[i], i);
		}

		// special cases: "" can be combine with any palindrome string
		if (map.containsKey("")) {
			int blankIdx = map.get("");
			for (int i = 0; i < words.length; i++) {
				if (isPalindrome(words[i])) {
					if (i == blankIdx)
						continue;
					res.add(Arrays.asList(blankIdx, i));
					res.add(Arrays.asList(i, blankIdx));
				}
			}
		}

		// find all string and reverse string pairs
		for (int i = 0; i < words.length; i++) {
			String cur_r = reverseStr(words[i]);
			if (map.containsKey(cur_r)) {
				int found = map.get(cur_r);
				if (found == i)
					continue;
				res.add(Arrays.asList(i, found));
			}
		}

		// find the pair s1, s2 that
		// case1 : s1[0:cut] is palindrome and s1[cut+1:] = reverse(s2) => (s2, s1)
		// case2 : s1[cut+1:] is palindrome and s1[0:cut] = reverse(s2) => (s1, s2)
		for (int i = 0; i < words.length; i++) {
			String cur = words[i];
			for (int cut = 1; cut < cur.length(); cut++) {
				if (isPalindrome(cur.substring(0, cut))) {
					String cut_r = reverseStr(cur.substring(cut));
					if (map.containsKey(cut_r)) {
						int found = map.get(cut_r);
						if (found == i)
							continue;
						res.add(Arrays.asList(found, i));
					}
				}
				if (isPalindrome(cur.substring(cut))) {
					String cut_r = reverseStr(cur.substring(0, cut));
					if (map.containsKey(cut_r)) {
						int found = map.get(cut_r);
						if (found == i)
							continue;
						res.add(Arrays.asList(i, found));
					}
				}
			}
		}
		// Collections.sort(res, Comparator.comparing(list -> list.get(0)));
		return res;
	}

	private boolean empty(int[] freq) {
		for (int f : freq)
			if (f > 0)
				return false;
		return true;
	}

	private String toString(int[] freq) {
		StringBuilder sb = new StringBuilder();
		char c = 'a';
		for (int f : freq) {
			while (f-- > 0)
				sb.append(c);
			c++;
		}
		return sb.toString();
	}

	public int minStickers(String[] stickers, String target) {
		// Optimization 1: Maintain frequency only for characters present in target
		int[] targetNaiveCount = new int[26];
		for (char c : target.toCharArray())
			targetNaiveCount[c - 'a']++;
		int[] index = new int[26];
		int N = 0; // no of distinct characters in target
		for (int i = 0; i < 26; i++)
			index[i] = targetNaiveCount[i] > 0 ? N++ : -1;
		int[] targetCount = new int[N];
		int t = 0;
		for (int c : targetNaiveCount)
			if (c > 0) {
				targetCount[t++] = c;
			}
		int[][] stickersCount = new int[stickers.length][N];
		for (int i = 0; i < stickers.length; i++) {
			for (char c : stickers[i].toCharArray()) {
				int j = index[c - 'a'];
				if (j >= 0)
					stickersCount[i][j]++;
			}
		}
		// Optimization 2: Remove stickers dominated by some other sticker
		int start = 0;
		for (int i = 0; i < stickers.length; i++) {
			for (int j = start; j < stickers.length; j++)
				if (j != i) {
					int k = 0;
					while (k < N && stickersCount[i][k] <= stickersCount[j][k])
						k++;
					if (k == N) {
						int[] tmp = stickersCount[i];
						stickersCount[i] = stickersCount[start];
						stickersCount[start++] = tmp;
						break;
					}
				}
		}
		// Perform BFS with target as source and an empty string as destination
		Queue<int[]> Q = new LinkedList<>();
		Set<String> visited = new HashSet<>();
		Q.add(targetCount);
		int steps = 0;
		while (!Q.isEmpty()) {
			steps++;
			int size = Q.size();
			while (size-- > 0) {
				int[] freq = Q.poll();
				String cur = toString(freq);
				if (visited.add(cur)) {
					// Optimization 3: Only use stickers that are capable of removing first
					// character from current string
					int first = cur.charAt(0) - 'a';
					for (int i = start; i < stickers.length; i++)
						if (stickersCount[i][first] != 0) {
							int[] next = freq.clone();
							for (int j = 0; j < N; j++)
								next[j] = Math.max(next[j] - stickersCount[i][j], 0);
							if (empty(next))
								return steps;
							Q.add(next);
						}
				}
			}
		}
		return -1;
	}

	public int kSimilarity(String s1, String tar) {
		Queue<String> q = new ArrayDeque<>();
		q.add(s1);

		int lvl = 0;
		while (q.size() > 0) {
			int size = q.size();
			while (size-- > 0) {
				String s = q.remove();
				if (s.equals(tar))
					return lvl;

				int i = 0;
				while (s.charAt(i) == tar.charAt(i))
					i++;

				int j = i;

				while (j < s.length()) {
					if (s.charAt(j) == tar.charAt(i) && tar.charAt(j) != s.charAt(j)) {
						StringBuilder sb = new StringBuilder(s);
						sb.setCharAt(i, s.charAt(j));
						sb.setCharAt(j, s.charAt(i));

						// A small optimization.
						if (sb.toString().equals(tar))
							return lvl + 1;

						q.add(sb.toString());
					}
					j++;
				}
			}
			lvl++;
		}
		return lvl;
	}

	// Driver code to test above
	public static void main(String args[]) {

	}
}

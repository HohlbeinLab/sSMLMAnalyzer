package com.wurgobes.sSMLMAnalyzer;

import java.util.Arrays;

public class levenshtein {

    public static int getTheClosestMatch(String[] strings, String target) {
        int distance = Integer.MAX_VALUE;
        int index = 0;
        for(int i = 0; i < strings.length; i ++){
            String compareString = strings[i];
            int currentDistance = levenshtein.calculate(compareString, target);
            if(currentDistance < distance) {
                distance = currentDistance;
                index = i;
            }
        }
        return index;
    }

    static int calculate(String x, String y) {
        if(y == null) y = "null";

        int[][] dp = new int[x.length() + 1][y.length() + 1];

        for (int i = 0; i <= x.length(); i++) {
            for (int j = 0; j <= y.length(); j++) {
                if (i == 0) {
                    dp[i][j] = j;
                }
                else if (j == 0) {
                    dp[i][j] = i;
                }
                else {
                    dp[i][j] = min(dp[i - 1][j - 1]
                                    + costOfSubstitution(x.charAt(i - 1), y.charAt(j - 1)),
                            dp[i - 1][j] + 1,
                            dp[i][j - 1] + 1);
                }
            }
        }

        return dp[x.length()][y.length()];
    }

    public static int costOfSubstitution(char a, char b) {
        return a == b ? 0 : 1;
    }

    public static int min(int... numbers) {
        return Arrays.stream(numbers)
                .min().orElse(Integer.MAX_VALUE);
    }
}
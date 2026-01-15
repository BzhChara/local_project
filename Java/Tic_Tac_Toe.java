import java.util.Arrays;
import java.util.Scanner;

public class Tic_Tac_Toe {
    public static int[] getIndex (String[][] board, String number) {
        for (int i = 0; i < board.length; i ++) {
            for (int j = 0; j < board[i].length; j ++) {
                if (board[i][j].equals(number)) {
                    return new int[] {i, j};
                }
            }
        }
        return new int[] {-1, -1};
    }

    public static void main(String[] args) {
        String[][] board = new String[][] {
            {"1", "2", "3"},
            {"4", "5", "6"},
            {"7","8", "9"}
        };

        //绘制游戏板
        System.out.println("The game board:");
        for (int i = 0; i < board.length; i++) {
            for (int j = 0; j < board[i].length; j++) {
                System.out.print(board[i][j] + " ");
                if ((j + 1) % 3 == 0) {
                    System.out.println();
                }
            }
        }

        System.out.println("The first one is 'O'!");
        System.out.println("The second one is 'X'!");

        Scanner sc = new Scanner(System.in);

        int time = 9;
        String winner = "nobody";

        //开始下棋
        while (true) {
            System.out.print("The first one's turn:");
            String first = sc.nextLine();
            while (true) {
                if (Arrays.equals(getIndex(board, first), new int[] {-1, -1})) {
                    System.out.print("Please input the correct number:");
                    first = sc.nextLine();
                }
                else {
                    board[getIndex(board, first)[0]][getIndex(board, first)[1]] = "O";
                    break;
                }
            }
            
            for (int i = 0; i < board.length; i++) {
                for (int j = 0; j < board[i].length; j++) {
                    System.out.print(board[i][j] + " ");
                    if ((j + 1) % 3 == 0) {
                        System.out.println();
                    }
                }
            }

            time = time - 1;
            
            //判断赢家
            for (int i =0; i < board.length; i++) {
                if (Arrays.equals(board[i], new String[] {"O", "O", "O"})) {
                    winner = "the first one";
                    break;
                }
                else if (Arrays.equals(board[i], new String[] {"X", "X", "X"})) {
                    winner = "the second one";
                    break;
                }
            }

            for (int j = 0; j < board[0].length; j++) {
                if (Arrays.equals(new String[] {board[0][j], board[1][j], board[2][j]}, new String[] {"O", "O", "O"})) {
                    winner = "the first one";
                    break;
                }
                else if (Arrays.equals(new String[] {board[0][j], board[1][j], board[2][j]}, new String[] {"X", "X", "X"})) {
                    winner = "the second one";
                    break;
                }
            }

            if (Arrays.equals(new String[] {board[0][0], board[1][1], board[2][2]}, new String[] {"O", "O", "O"})) {
                winner = "the first one";
                break;
            }
            else if (Arrays.equals(new String[] {board[0][2], board[1][1], board[2][0]}, new String[] {"O", "O", "O"})) {
                winner = "the first one";
                break;
            }
            else if (Arrays.equals(new String[] {board[0][0], board[1][1], board[2][2]}, new String[] {"X", "X", "X"})) {
                winner = "the second one";
                break;
            }
            else if (Arrays.equals(new String[] {board[0][2], board[1][1], board[2][0]}, new String[] {"X", "X", "X"})) {
                winner = "the second one";
                break;
            }

            if (!winner.equals("nobody")) {
                break;
            }

            //游戏轮数耗尽
            if (time <= 0) {
                System.out.println("Draw, game over!");
                break;
            }
            
            System.out.print("The second one's turn:");
            String second = sc.nextLine();
            while (true) {
                if (Arrays.equals(getIndex(board, second), new int[] {-1, -1})) {
                    System.out.print("Please input the correct number:");
                    second = sc.nextLine();
                }
                else {
                    board[getIndex(board, second)[0]][getIndex(board, second)[1]] = "X";
                    break;
                }
            }  

            for (int i = 0; i < board.length; i++) {
                for (int j = 0; j < board[i].length; j++) {
                    System.out.print(board[i][j] + " ");
                    if ((j + 1) % 3 == 0) {
                        System.out.println();
                    }
                }
            }

            time = time - 1;

            //判断赢家
            for (int i =0; i < board.length; i++) {
                if (Arrays.equals(board[i], new String[] {"O", "O", "O"})) {
                    winner = "the first one";
                    break;
                }
                else if (Arrays.equals(board[i], new String[] {"X", "X", "X"})) {
                    winner = "the second one";
                    break;
                }
            }

            for (int j = 0; j < board[0].length; j++) {
                if (Arrays.equals(new String[] {board[0][j], board[1][j], board[2][j]}, new String[] {"O", "O", "O"})) {
                    winner = "the first one";
                    break;
                }
                else if (Arrays.equals(new String[] {board[0][j], board[1][j], board[2][j]}, new String[] {"X", "X", "X"})) {
                    winner = "the second one";
                    break;
                }
            }

            if (Arrays.equals(new String[] {board[0][0], board[1][1], board[2][2]}, new String[] {"O", "O", "O"})) {
                winner = "the first one";
                break;
            }
            else if (Arrays.equals(new String[] {board[0][2], board[1][1], board[2][0]}, new String[] {"O", "O", "O"})) {
                winner = "the first one";
                break;
            }
            else if (Arrays.equals(new String[] {board[0][0], board[1][1], board[2][2]}, new String[] {"X", "X", "X"})) {
                winner = "the second one";
                break;
            }
            else if (Arrays.equals(new String[] {board[0][2], board[1][1], board[2][0]}, new String[] {"X", "X", "X"})) {
                winner = "the second one";
                break;
            }

            if (!winner.equals("nobody")) {
                break;
            }

            //游戏轮数耗尽
            if (time <= 0) {
                System.out.println("Draw, game over!");
                break;
            }
        }
        System.out.println("The winner is " + winner + " !");
        sc.close();
    }
}
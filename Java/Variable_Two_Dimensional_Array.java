import java.util.Arrays;
import java.util.Scanner;

public class Variable_Two_Dimensional_Array {
    public static void main(String[] args){
      
        Scanner sc = new Scanner(System.in);

        System.out.println("请输入要处理的数据组数: ");
        int epoch = sc.nextInt();
        sc.nextLine();
        int[][] array = new int[epoch][];

        for (int i = 0; i < epoch; i++){
            System.out.print("请输入第 " + (i+1) + " 组数据（空格分隔）: ");
            String str = sc.nextLine();

            String[] strArray = str.split(" ");
            int[] tmp = new int[strArray.length];

            for (int j = 0; j < strArray.length; j++){
                tmp[j] = Integer.parseInt(strArray[j]);
            }

            array[i] = tmp;
        }

        System.out.println("处理后的二维数组:");
        for (int i = 0; i < array.length; i++) {
            System.out.println("第 " + i + " 行: " + Arrays.toString(array[i]));
        }

        sc.close();
    }
}
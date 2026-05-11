package Main;

import java.util.List;

public class Main {
    public static void main(String[] args) {
        List<Integer> list = List.of(12, 34, 56);
        Integer[] array = list.toArray(new Integer[0]);
        for (Integer n : array) {
            System.out.println(n);
        }
    }
}

package HashMap_TreeMap_Implement;

import java.util.*;

public class HashMap_TreeMap_Implement {
    public static void main(String[] args) {
        List<Student> students = List.of(
                new Student("Tom", 77),
                new Student("Bob", 66),
                new Student("Lily", 99),
                new Student("Alice", 85),
                new Student("Jack", 59),
                new Student("Rose", 92),
                new Student("Mike", 70));

        // TODO 1:
        // 创建一个 HashMap<Student, Integer>，变量名叫 scoreMap
        // key 是 Student，value 是 score
        Map<Student, Integer> scoreMap = new HashMap<>();

        // TODO 2:
        // 遍历 students，把每个学生和他的分数放入 scoreMap
        // 提示：scoreMap.put(student, student.score);
        for (Student student : students) {
            scoreMap.put(student, student.score);
        }

        // TODO 3:
        // 测试 HashMap 查询：
        // 用 new Student("Bob", 66) 查询分数，期望输出 66
        Integer bobScore = scoreMap.get(new Student("Bob", 66));
        System.out.println("Bob score = " + bobScore);

        // TODO 4:
        // 创建一个 TreeMap<String, Integer>，变量名叫 rangeMap
        // 要求分数段从高到低排序：
        // 90-100, 80-89, 70-79, 60-69, 0-59
        //
        // 提示：
        // new TreeMap<>(new Comparator<String>() { ... })
        Map<String, Integer> rangeMap = new TreeMap<>(new Comparator<String>() {
            @Override
            public int compare(String s1, String s2) {
                // 90-100 > 80-89 > 70-79 > 60-69 > 0-59
                // 可以通过比较字符串的第一个字符来实现排序
                char c1 = s1.charAt(0);
                char c2 = s2.charAt(0);
                return Character.compare(c2, c1); // 降序排序
            }
        });

        // TODO 5:
        // 遍历 students，根据每个学生的 score 得到分数段
        // 然后统计每个分数段的人数
        //
        // 提示：
        // String range = getScoreRange(student.score);
        // Integer count = rangeMap.get(range);
        // 如果 count == null，说明这个分数段还没人，放入 1
        // 否则放入 count + 1
        for (Student student : students) {
            String range = getScoreRange(student.score);
            Integer count = rangeMap.get(range);
            if (count == null) {
                rangeMap.put(range, 1);
            } else {
                rangeMap.put(range, count + 1);
            }
        }

        // 输出统计结果
        System.out.println("Score ranges:");
        for (Map.Entry<String, Integer> entry : rangeMap.entrySet()) {
            System.out.println(entry.getKey() + " = " + entry.getValue());
        }
    }

    static String getScoreRange(int score) {
        if (score >= 90) {
            return "90-100";
        }
        if (score >= 80) {
            return "80-89";
        }
        if (score >= 70) {
            return "70-79";
        }
        if (score >= 60) {
            return "60-69";
        }
        return "0-59";
    }
}

class Student {
    String name;
    int score;

    Student(String name, int score) {
        this.name = name;
        this.score = score;
    }

    @Override
    public String toString() {
        return "{Student: name=" + name + ", score=" + score + "}";
    }

    // TODO 6:
    // 重写 equals
    // 判断规则：name 相同，并且 score 相同，就认为是同一个 Student
    //
    // 提示：
    // if (this == o) return true;
    // if (o instanceof Student s) { ... }
    // return false;
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o instanceof Student s) {
            return Objects.equals(this.name, s.name) && this.score == s.score;
        }
        return false;
    }

    // TODO 7:
    // 重写 hashCode
    // 要和 equals 保持一致
    //
    // 提示：
    // return Objects.hash(name, score);
    @Override
    public int hashCode() {
        return Objects.hash(name, score);
    }
}
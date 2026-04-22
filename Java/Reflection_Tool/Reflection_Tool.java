package Reflection_Tool;

import java.lang.reflect.Field;

public class Reflection_Tool {
    public static void main(String[] args) throws Exception {
        // 1. 准备一个具体的对象实例（小明）
        Person p = new Person("Xiao Ming", 20);

        // 2. 调用我们的通用打印工具
        printAllFields(p);
    }

    /**
     * 通用工具：能够透视任何对象的私有字段
     */
    public static void printAllFields(Object obj) throws Exception {
        // 获取该对象的 Class 实例（获取“说明书”）
        Class<?> cls = obj.getClass();

        // 获取该类自己声明的所有字段（包括 private）
        Field[] fields = cls.getDeclaredFields();

        System.out.println("开始透视类: " + cls.getName());

        for (Field f : fields) {
            // 【关键步骤】强行开启私有访问权限（打破 private 限制）
            f.setAccessible(true);

            // 获取字段名
            String name = f.getName();

            // 获取该字段在具体对象 obj 中的值
            Object value = f.get(obj);

            System.out.println("字段名: " + name + " | 对应的值: " + value);
        }
    }
}

class Person {
    @SuppressWarnings("unused") // 专门抑制“未使用”警告
    private String name;

    @SuppressWarnings("unused")
    private int age;

    public Person(String name, int age) {
        this.name = name;
        this.age = age;
    }
}
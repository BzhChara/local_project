package PriorityQueue_Exercises;

import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.Queue;

public class PriorityQueue_Exercises {
    public static void main(String[] args) {

        Queue<User> queue = new PriorityQueue<>(new Comparator<User>() {
            @Override
            public int compare(User u1, User u2) {
                // TODO 1:
                // VIP 等级高的排前面
                // 提示：
                // u2.vipLevel 和 u1.vipLevel 比较，因为 PriorityQueue 默认“小的”在前面
                // 如果想让 vipLevel 大的先出队，要反过来比较
                if (u1.vipLevel != u2.vipLevel) {
                    return Integer.compare(u2.vipLevel, u1.vipLevel);
                }

                return Integer.compare(u1.arrivalOrder, u2.arrivalOrder);
            }
        });

        queue.offer(new User("普通用户A", 0, 1));
        queue.offer(new User("VIP1用户B", 1, 2));
        queue.offer(new User("普通用户C", 0, 3));
        queue.offer(new User("VIP3用户D", 3, 4));
        queue.offer(new User("VIP2用户E", 2, 5));
        queue.offer(new User("VIP3用户F", 3, 6));
        queue.offer(new User("VIP1用户G", 1, 7));

        System.out.println("开始服务：");

        // TODO 3:
        // 使用 while 循环，不断从 queue 中取出用户
        // 提示：
        // queue.isEmpty()
        // queue.poll()
        while (!queue.isEmpty()) {
            System.out.println(queue.poll());
        }
    }
}

class User {
    String name;
    int vipLevel;
    int arrivalOrder;

    User(String name, int vipLevel, int arrivalOrder) {
        this.name = name;
        this.vipLevel = vipLevel;
        this.arrivalOrder = arrivalOrder;
    }

    @Override
    public String toString() {
        return name + "，VIP等级=" + vipLevel + "，到达顺序=" + arrivalOrder;
    }
}
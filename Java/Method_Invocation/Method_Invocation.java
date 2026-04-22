package Method_Invocation;

import java.lang.reflect.Method;

public class Method_Invocation {

	public static void main(String[] args) throws Exception {
		String name = "Xiao Ming";
		int age = 20;
		Person p = new Person();
		Method[] methods = p.getClass().getMethods();
		for (Method m : methods) {
			if ("setName".equals(m.getName())) {
				m.invoke(p, name);
			} else if ("setAge".equals(m.getName())) {
				m.invoke(p, age);
			} else {
				continue;
			}
		}

		System.out.println(p.getName()); // "Xiao Ming"
		System.out.println(p.getAge()); // 20
	}
}

class Person {
	private String name;
	private int age;

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public int getAge() {
		return age;
	}

	public void setAge(int age) {
		this.age = age;
	}
}
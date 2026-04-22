package Annotation_Reflection;

import java.lang.annotation.*;
import java.lang.reflect.*;

// 1. 定义注解：标记哪些参数不能为负数
@Target(ElementType.PARAMETER)
@Retention(RetentionPolicy.RUNTIME)
@interface NotNegative {
}

// 2. 定义业务接口：模拟你的校准逻辑
interface Calibrator {
    void calibrate(String sensorName, @NotNegative int value);
}

// 3. 编写真实的业务实现类（只管写业务，不管校验）
class RealCalibrator implements Calibrator {
    @Override
    public void calibrate(String sensorName, int value) {
        System.out.println(">>> [业务执行] 传感器 " + sensorName + " 校准成功，数值为: " + value);
    }
}

// 4. 编写代理处理器：这是“管家”，负责拦截并检查标签
class ValidationHandler implements InvocationHandler {
    private final Object realObj;

    public ValidationHandler(Object realObj) {
        this.realObj = realObj;
    }

    @Override
    public Object invoke(Object proxy, Method method, Object[] args) throws Throwable {
        // 获取该方法所有参数的注解（二维数组）
        Annotation[][] parameterAnnotations = method.getParameterAnnotations();

        for (int i = 0; i < parameterAnnotations.length; i++) {
            for (Annotation anno : parameterAnnotations[i]) {
                // 如果发现 @NotNegative 标签
                if (anno instanceof NotNegative) {
                    int val = (int) args[i]; // 获取对应的参数值
                    if (val < 0) {
                        // 发现违规，直接抛出异常，不再执行后面的 method.invoke
                        System.err.println("!!! [拦截成功] 发现非法参数: " + val);
                        throw new IllegalArgumentException("报错：校准值不能为负数！");
                    }
                }
            }
        }

        // 校验通过，放行，执行真实的业务代码
        return method.invoke(realObj, args);
    }
}

// 5. 测试类
public class Annotation_Reflection {
    public static void main(String[] args) {
        // 创建真实对象
        Calibrator realWorker = new RealCalibrator();

        // 创建代理对象
        Calibrator proxy = (Calibrator) Proxy.newProxyInstance(
                Calibrator.class.getClassLoader(),
                new Class[] { Calibrator.class },
                new ValidationHandler(realWorker));

        // --- 测试场景 1：正常数据 ---
        System.out.println("--- 尝试输入正常数据 ---");
        proxy.calibrate("氮肥传感器", 100);

        // --- 测试场景 2：非法数据 ---
        System.out.println("\n--- 尝试输入非法数据 ---");
        try {
            proxy.calibrate("磷肥传感器", -5);
        } catch (Exception e) {
            System.out.println("捕获到异常: " + e.getMessage());
        }
    }
}
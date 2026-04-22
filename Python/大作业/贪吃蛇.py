import os
import sys
import time
import random
import pygame

#如果缺少pygame请使用pip install pygame安装

#每次游玩会记录你的最好成绩

#初始化pygame
pygame.init()

#游戏窗口和基本参数
screen_width, screen_height = input("请输入分辨率(例如1920*1080):").split('*')
screen_width = int(screen_width)
screen_height = int(screen_height)

#定义字体和字号
font = pygame.font.SysFont("Arial Black", 36)

#蛇的初始大小和速度
snake_size = 20
snake_speed = 10

#定义颜色
white = (255, 255, 255)
red = (255, 0, 0)

#初始化窗口
screen = pygame.display.set_mode((screen_width, screen_height))
pygame.display.set_caption("Snake Game")#标题(全屏就看不到了。。。)

#初始蛇身和方向
snake = [(100, 100), (90, 100), (80, 100)]
snake_direction = (1, 0)

#初始食物位置
food = (random.randint(1, (screen_width - snake_size) // snake_size) * snake_size,#计算能够容纳整数个蛇身块的最大宽度
        random.randint(1, (screen_height - snake_size) // snake_size) * snake_size)

#定义分数
score = 0
if  not (os.path.exists("./score.txt")):
    score_path = open("./score.txt", 'w')
    score_path.write(str(score))
    score_path.close()

#游戏BGM
file_path = "TheFatRat - Fire.flac"#打开文件(可自行替换文件)
pygame.mixer.init()
pygame.mixer.music.load(file_path)#加载文件
pygame.mixer.music.play()#播放文件

#游戏主循环
while True:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit()

    #响应输入(上下左右)
    keys = pygame.key.get_pressed()
    for key in keys:
        if (keys[pygame.K_UP] and snake_direction != (0, 1)):
            snake_direction = (0, -1)
        elif (keys[pygame.K_DOWN] and snake_direction != (0, -1)):
            snake_direction = (0, 1)
        elif (keys[pygame.K_LEFT] and snake_direction != (1, 0)):
            snake_direction = (-1, 0)
        elif (keys[pygame.K_RIGHT] and snake_direction != (-1, 0)):
            snake_direction = (1, 0)

    #更新蛇的位置
    snake_head = (snake[0][0] + snake_direction[0] * snake_size, snake[0][1] + snake_direction[1] * snake_size)
    snake.insert(0, snake_head)

    #判断是否吃到食物
    if (snake_head == food):
        #每吃到一次食物速度加快(设置上限)
        #想要添加挑战性可自行删除
        if (snake_speed <= 30):
            snake_speed = snake_speed + 5
        score = score + 1
        #重新随机生成食物
        food = (random.randint(1, (screen_width - snake_size) // snake_size) * snake_size,#计算能够容纳整数个蛇身块的最大宽度
                random.randint(1, (screen_height - snake_size) // snake_size) * snake_size)
    #如果没有吃到食物，则删除蛇尾，实现移动效果
    else:
        snake.pop()

    #判断是否碰到墙壁或自身
    if (snake_head[0] < 0 or snake_head[0] >= screen_width or snake_head[1] < 0 or snake_head[1] >= screen_height or snake_head in snake[1:]):
        #读取分数
        score_path = open("./score.txt", 'r+')
        tmp = score_path.readline()
        #print(tmp)
        if score > int(tmp):
            #文件定位
            score_path.seek(0)
            #清空文件
            score_path.truncate()
            #写入分数
            score_path.write(str(score))
            score_path.close()
            #渲染文字
            text_1 = font.render("Game Over!", True, white)
            text_2 = font.render("New Score:"+str(score), True, white)
        else:
            #渲染文字
            text_1 = font.render("Game Over!", True, white)
            text_2 = font.render("Score:"+str(score), True, white)
            score_path.close()
        #清空屏幕
        screen.fill((0, 0, 0))
        #将文字绘制到窗口
        screen.blit(text_1, (screen_width // 2.3, screen_height // 2.1))
        screen.blit(text_2, (screen_width // 2.3, screen_height // 1.9))
        #更新屏幕
        pygame.display.update()
        #淡出音乐
        pygame.mixer.music.fadeout(8000) 
        time.sleep(8)
        
        pygame.quit()
        sys.exit()

    # 检查音乐是否结束
    if not pygame.mixer.music.get_busy():
        #读取分数
        score_path = open("./score.txt", 'r+')
        tmp = score_path.readline()
        #print(tmp)
        if score > int(tmp):
            #文件定位
            score_path.seek(0)
            #清空文件
            score_path.truncate()
            #写入分数
            score_path.write(str(score))
            score_path.close()
            #渲染文字
            text_1 = font.render("Game Over!", True, white)
            text_2 = font.render("New Score:"+str(score), True, white)
        else:
            #渲染文字
            text_1 = font.render("Game Over!", True, white)
            text_2 = font.render("Score:"+str(score), True, white)
            score_path.close()
        #清空屏幕
        screen.fill((0, 0, 0))
        #将文字绘制到窗口
        screen.blit(text_1, (screen_width // 2.3, screen_height // 2.1))
        screen.blit(text_2, (screen_width // 2.2, screen_height // 1.9))
        #更新屏幕
        pygame.display.update()
        #淡出音乐
        pygame.mixer.music.fadeout(8000) 
        time.sleep(8)
        
        pygame.quit()
        sys.exit()

    #绘制蛇和食物
    screen.fill((0, 0, 0))
    pygame.draw.rect(screen, red, (food[0], food[1], snake_size, snake_size))
    for segment in snake:
        pygame.draw.rect(screen, white, (segment[0], segment[1], snake_size, snake_size))
    pygame.display.flip()
    #设置游戏速度
    pygame.time.Clock().tick(snake_speed)
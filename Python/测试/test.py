comment = """I shall be telling this with a sigh
    Somewhere ages and ages hence:
    Two roads diverged in a wood, and I--
    I took the one less traveled by,
    And that has made all the difference."""
tmp = comment.split("\n")
for line in tmp:
    print(line.lstrip())

string = "Hello hello!"
print(len(string)) # 空格算字符串

num1 = 4*3**2 +1
num2 = -4**2
print(num1, num2)

arr1 = input().split(" ")
print(arr1, type(arr1))
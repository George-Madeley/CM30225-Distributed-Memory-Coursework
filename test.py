
num_list = []

for i in range(45):
  if i not in num_list:
    num_list.append(i)
  if 2 * i not in num_list:
    num_list.append(2 * i)
  if 3 * i not in num_list:
    num_list.append(3 * i)
  if 4 * i not in num_list:
    num_list.append(4 * i)
num_list.sort()
print(num_list)
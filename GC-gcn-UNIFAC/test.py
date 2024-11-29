import random

# 원본 리스트
lst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# 인덱스를 섞기
indices = list(range(len(lst)))
random.shuffle(indices)

# 섞인 인덱스로 리스트 요소 호출
shuffled_lst = [lst[i] for i in indices]

print("원본 리스트:", lst)
print("섞인 인덱스:", indices)
print("섞인 순서의 리스트:", shuffled_lst)

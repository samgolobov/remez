# remez


if my < x1[0]:
    if (sign(f(my) - p(my)) == sign(f(x1[0]) - p(x1[0]))):
        x1[0] = my
    else:
        x1[-1] = my
elif my > x1[-1]:
    if (sign(f(my) - p(my)) == sign(f(x1[0]) - p(x1[0]))):
        x1[-1] = my
    else:
        x1[0] = my
else:
    if sign(f(my) - p(my)) == sign(f(x1[my_Index - 1]) - p(x1[x1[my_Index - 1]])):
        x1[my_Index - 1] = my
    else:
        x1[my_Index + 1] = my

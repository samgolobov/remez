if my<x[0]:
    if sign(testf(my)-p(my)) == sign(testf(x[0])-p(x[0])):
        x[0]=my
    else: 
        x[-1]=my
if my>x[-1]:
    if sign(testf(my)-p(my)) == sign(testf(x[-1])-p(x[-1])):
        x[-1]=my
    else:
        x[0]=my
else:
    for i in range(len(x)):
        if x[i]>=my:
            x_i = i-1
            break
    if sign(testf(my)-p(my)) == sign(testf(x[x_i])-p(x[x_i])):
        x[x_i] = my
    else:
        x[x_i + 1] = my


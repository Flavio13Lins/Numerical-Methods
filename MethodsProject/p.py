from sympy.parsing.sympy_parser import parse_expr


def aplica_fyt(y, t, fn):
    return fn.subs('y', y).subs('t', t)


def calc_eul(y, t, h, fn):
    return y + h * aplica_fyt(y, t, fn)


def euler_retorno_n_primeiros(param):
    [y0, t0, h, _, func, ordem] = param.split()[1:]
    ordem = int(ordem)
    h = float(h)
    y0 = float(y0)
    t0 = float(t0)
    fn = parse_expr(func)
    valores = list(range(ordem))
    for x in range(ordem):
        y0 = calc_eul(y0, t0, h, fn)
        t0 = t0 + h
        valores[x] = y0
    return valores


def euler(entr_euler):
    [y0, t0, pas, func] = entr_euler.split()[1:]
    print('Metodo de Euler')
    h = 0.1
    y0 = float(y0)
    t0 = float(t0)
    pas = int(pas)
    fn = parse_expr(func)
    print('y( {} ) = {}'.format(t0, y0))
    print('h = {}'.format(h))
    for x in range(0, pas+1):
        print('{} {}'.format(x, y0))
        y0 = calc_eul(y0, t0, h, fn)
        t0 = t0 + h
    return


'''
def calc_euler_inv(y, t, h, fn):
    return


def euler_inverso_retorno_n_primeiros(param):
    # adam_bashforth_by_euler_inverso 0 0 0.1 20 1-t+4*y 6
    return
'''


def euler_inverso(entr_euler_inv):
    [y0, t0, h, pas, func] = entr_euler_inv.split()[1:]
    print('Metodo de Euler Inverso')
    y0 = float(y0)
    t0 = float(t0)
    h = float(h)
    pas = int(pas)
    fn = parse_expr(func)
    print('y( {} ) = {}'.format(t0, y0))
    print('h = {}'.format(h))
    for x in range(0, pas+1):
        print('{} {}'.format(x, y0))
        # y0 = sympy.solvers.solve(y0 + (h * fn.subs(t, t0+h)) - 'y', 'y')[0]  # aproximaçao melhor como do livro
        y0 = y0 + h * aplica_fyt(calc_eul(y0, t0, h, fn), t0 + h, fn)  # aproximaçao como pdf dos monitores
        t0 = t0 + h
    return


def euler_aprimorado(entr_euler_apri):
    [y0, t0, h, pas, func] = entr_euler_apri.split()[1:]
    print('Metodo de Euler Aprimorado')
    y0 = float(y0)
    t0 = float(t0)
    h = float(h)
    pas = int(pas)
    k = parse_expr(func)
    print('y( {} ) = {}'.format(t0, y0))
    print('h = {}'.format(h))
    for x in range(0, pas+1):
        print('{} {}'.format(x, y0))
        k1 = aplica_fyt(y0, t0, k)
        k2 = aplica_fyt(y0 + h * k1, t0 + h, k)
        y0 = y0 + (h/2) * (k1 + k2)
        t0 = t0 + h
    return


def runge_kutta(entr_runge_kutta):
    [y0, t0, h, pas, func] = entr_runge_kutta.split()[1:]
    print('Metodo de Runge-Kutta')
    y0 = float(y0)
    t0 = float(t0)
    h = float(h)
    pas = int(pas)
    k = parse_expr(func)
    print('y( {} ) = {}'.format(t0, y0))
    print('h = {}'.format(h))
    for x in range(0, pas+1):
        print('{} {}'.format(x, y0))
        k1 = aplica_fyt(y0, t0, k)
        k2 = aplica_fyt(y0 + h * k1/2, t0 + h/2, k)
        k3 = aplica_fyt(y0 + h * k2/2, t0 + h/2, k)
        k4 = aplica_fyt(y0 + h * k3, t0 + h, k)
        y0 = y0 + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
        t0 = t0 + h
    return


def aplica_bf2_fyt(y0, y1, t0, t1, fn):
    return (3 / 2) * aplica_fyt(y1, t1, fn) \
         - (1 / 2) * aplica_fyt(y0, t0, fn)


def aplica_bf3_fyt(y0, y1, y2, t0, t1, t2, fn):
    return (23 / 12) * aplica_fyt(y2, t2, fn) \
            - (4 / 3) * aplica_fyt(y1, t1, fn) \
            + (5 / 12) * aplica_fyt(y0, t0, fn)


def aplica_bf4_fyt(y0, y1, y2, y3, t0, t1, t2, t3, fn):
    return (55 / 24) * aplica_fyt(y3, t3, fn) \
            - (59 / 24) * aplica_fyt(y2, t2, fn) \
            + (37 / 24) * aplica_fyt(y1, t1, fn) \
            - (3 / 8) * aplica_fyt(y0, t0, fn)


def aplica_bf5_fyt(y0, y1, y2, y3, y4, t0, t1, t2, t3, t4, fn):
    return (1901 / 720) * aplica_fyt(y4, t4, fn) \
            - (1387 / 360) * aplica_fyt(y3, t3, fn) \
            + (109 / 30) * aplica_fyt(y2, t2, fn) \
            - (637 / 360) * aplica_fyt(y1, t1, fn) \
            + (251 / 720) * aplica_fyt(y0, t0, fn)


def aplica_bf6_fyt(y0, y1, y2, y3, y4, y5, t0, t1, t2, t3, t4, t5, fn):
    return (4277 / 1440) * aplica_fyt(y5, t5, fn) \
            - (2641 / 480) * aplica_fyt(y4, t4, fn) \
            + (4991 / 720) * aplica_fyt(y3, t3, fn) \
            - (3649 / 720) * aplica_fyt(y2, t2, fn) \
            + (959 / 480) * aplica_fyt(y1, t1, fn) \
            - (95 / 288) * aplica_fyt(y0, t0, fn)


def aplica_bf7_fyt(y0, y1, y2, y3, y4, y5, y6, t0, t1, t2, t3, t4, t5, t6, fn):
    return (198721 / 60480) * aplica_fyt(y6, t6, fn) \
            - (18637 / 2520) * aplica_fyt(y5, t5, fn) \
            + (235183 / 20160) * aplica_fyt(y4, t4, fn) \
            - (10754 / 945) * aplica_fyt(y3, t3, fn) \
            + (135713 / 20160) * aplica_fyt(y2, t2, fn) \
            - (5603 / 2520) * aplica_fyt(y1, t1, fn) \
            + (19087 / 60480) * aplica_fyt(y0, t0, fn)


def aplica_bf8_fyt(y0, y1, y2, y3, y4, y5, y6, y7, t0, t1, t2, t3, t4, t5, t6, t7, fn):
    return (16083 / 4480) * aplica_fyt(y7, t7, fn) \
            - (1152169 / 120960) * aplica_fyt(y6, t6, fn) \
            + (242653 / 13440) * aplica_fyt(y5, t5, fn) \
            - (296053 / 13440) * aplica_fyt(y4, t4, fn) \
            + (2102243 / 120960) * aplica_fyt(y3, t3, fn) \
            - (115747 / 13440) * aplica_fyt(y2, t2, fn) \
            + (32863 / 13440) * aplica_fyt(y1, t1, fn) \
            - (5257 / 17280) * aplica_fyt(y0, t0, fn)


def prep_param_to_bashforth(met, y0, y, t, h, pas, func, ordem):
    ordem = int(ordem)
    ys = list(range(ordem+1))
    ys[0] = str(y0)
    for x in range(ordem):
        ys[x+1] = str(y[x])
    t = str(t)
    h = str(h)
    pas = str(pas)
    func = str(func)
    param_input = str(met)
    for x in range(ordem+1):
        param_input = param_input + ' ' + ys[x]
    ordem = str(ordem)
    param_input = param_input + ' ' + t + ' ' + h + ' ' + pas + ' ' + func + ' ' + ordem
    return param_input


def adam_bashforth_execute_out(param):
    [t0, h, pas, func, ordem] = param.split()[-5:]
    t0 = float(t0)
    h = float(h)
    pas = int(pas)
    k = parse_expr(func)
    ordem = int(ordem)
    y = list(range(pas + 1))
    for x in range(ordem):
        y[x] = float(param.split()[1 + x])
    print('y( {} ) = {}'.format(t0, y[0]))
    print('h = {}'.format(h))
    t0 = t0 + (ordem - 1) * h
    if ordem == 1:
        print('0 {}'.format(y[0]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_fyt(y[x], t0, k)
            t0 = t0 + h
            print('{} {}'.format(x, y[x]))
    if ordem == 2:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_bf2_fyt(y[x - 1], y[x - 2],
                                                 t0, t0 + h, k)
            t0 = t0 + h
            print('{} {}'.format(x, y[x]))
    if ordem == 3:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_bf3_fyt(y[x], y[x - 1], y[x - 2],
                                                 t0, t0 + h, t0 + 2 * h, k)
            t0 = t0 + h
            print('{} {}'.format(x, y[x]))
    if ordem == 4:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_bf4_fyt(y[x], y[x - 1], y[x - 2], y[x - 3],
                                                 t0, t0 + h, t0 + 2 * h, t0 + 3 * h, k)
            t0 = t0 + h
            print('{} {}'.format(x, y[x]))
    if ordem == 5:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_bf5_fyt(y[x], y[x - 1], y[x - 2], y[x - 3], y[x - 4],
                                                 t0, t0 + h, t0 + 2 * h, t0 + 3 * h, t0 + 4 * h, k)
            t0 = t0 + h
            print('{} {}'.format(x, y[x]))
    if ordem == 6:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_bf6_fyt(y[x], y[x - 1], y[x - 2], y[x - 3], y[x - 4], y[x - 5],
                                                 t0, t0 + h, t0 + 2 * h, t0 + 3 * h, t0 + 4 * h, t0 + 5 * h, k)
            t0 = t0 + h
            print('{} {}'.format(x, y[x]))
    if ordem == 7:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_bf7_fyt(y[x], y[x - 1], y[x - 2], y[x - 3], y[x - 4], y[x - 5], y[x - 6],
                                                 t0, t0 + h, t0 + 2 * h, t0 + 3 * h, t0 + 4 * h, t0 + 5 * h, t0 + 6 * h,
                                                 k)
            t0 = t0 + h
            print('{} {}'.format(x, y[x]))
    if ordem == 8:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_bf8_fyt(y[x], y[x - 1], y[x - 2], y[x - 3], y[x - 4], y[x - 5], y[x - 6],
                                                 y[x - 7],
                                                 t0, t0 + h, t0 + 2 * h, t0 + 3 * h, t0 + 4 * h, t0 + 5 * h, t0 + 6 * h,
                                                 t0 + 7 * h, k)
            t0 = t0 + h
            print('{} {}'.format(x, y[x]))
    return


def adam_bashforth(entrada_adam_bashforth):
    print('Metodo Adan-Bashforth')
    adam_bashforth_execute_out(entrada_adam_bashforth)
    return


def adam_bashforth_by_euler(entrada_adam_bashforth_euler):
    print('Metodo Adan-Bashforth por Euler')
    [met, y0, t0, h, pas, func, ordem] = entrada_adam_bashforth_euler.split()
    ordem = int(ordem)
    y1 = euler_retorno_n_primeiros(entrada_adam_bashforth_euler)
    orded_param = prep_param_to_bashforth(met, y0, y1, t0, h, pas, func, ordem)
    adam_bashforth_execute_out(orded_param)
    return


'''
def adam_bashforth_by_euler_inverso(entrada_adam_bashforth_euler_inverso):
    print('Metodo Adan-Bashforth por Euler Inverso')
    [metodo, y0, t0, h, pas, func, ordem] = entrada_adam_bashforth_euler_inverso.split()
    return
'''


entrada = str(input())
metodo = entrada.split()[0]
if metodo == 'euler':
    euler(entrada)
elif metodo == 'euler_inverso':
    euler_inverso(entrada)
elif metodo == 'euler_aprimorado':
    euler_aprimorado(entrada)
elif metodo == 'runge_kutta':
    runge_kutta(entrada)
elif metodo == 'adam_bashforth':
    adam_bashforth(entrada)
elif metodo == 'adam_bashforth_by_euler':
    adam_bashforth_by_euler(entrada)
elif metodo == 'adam_bashforth_by_euler_inverso':
    # todo adam_bashforth_by_euler_inverso(entrada)
elif metodo == 'adam_bashforth_by_euler_aprimorado':
    # todo adam_bashforth_by_euler_aprimorado(entrada)
elif metodo == 'adam_bashforth_by_runge_kutta':
    # todo adam_bashforth_by_runge_kutta(entrada)
elif metodo == 'adam_multon':
    # todo adam_multon(entrada)
elif metodo == 'adam_multon_by_euler':
    # todo adam_multon_by_euler(entrada)
elif metodo == 'adam_multon_by_euler_inverso':
    # todo adam_multon_by_euler_inverso(entrada)
elif metodo == 'adam_multon_by_euler_aprimorado':
    # todo adam_multon_by_euler_aprimorado(entrada)
elif metodo == 'adam_multon_by_runge_kutta':
    # todo adam_multon_by_runge_kutta(entrada)
else:
    print('.')


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
    valseul = list(range(ordem))
    for x in range(ordem):
        y0 = calc_eul(y0, t0, h, fn)
        t0 = t0 + h
        valseul[x] = y0
    return valseul


def euler(entr_euler):
    [y0, t0, pas, func] = entr_euler.split()[1:]
    print('Metodo de Euler')
    h = 0.1
    y0 = float(y0)
    t0 = float(t0)
    pas = int(pas)
    fn = parse_expr(func)
    print('y( {} ) = {}\nh = {}'.format(t0, y0, h))
    for x in range(0, pas+1):
        print('{} {}'.format(x, y0))
        y0 = calc_eul(y0, t0, h, fn)
        t0 = t0 + h
    return


def calc_eul_inv(y, t, h, fn):
    print('inv', y, t, h, fn, 'inv')
    print('r', y + h * aplica_fyt(calc_eul(y, t, h, fn), t + h, fn))
    return y + h * aplica_fyt(calc_eul(y, t, h, fn), t + h, fn)


def euler_inverso_retorno_n_primeiros(param):
    # adam_bashforth_by_euler_inverso 0 0 0.1 20 1-t+4*y 6
    [y0, t0, h, _, func, ordem] = param.split()[1:]
    ordem = int(ordem)
    h = float(h)
    y0 = float(y0)
    t0 = float(t0)
    fn = parse_expr(func)
    valsinv = list(range(ordem))
    for x in range(ordem):
        y0 = calc_eul_inv(y0, t0, h, fn)
        t0 = t0 + h
        valsinv[x] = y0
    return valsinv


def euler_inverso(entr_euler_inv):
    [y0, t0, h, pas, func] = entr_euler_inv.split()[1:]
    print('Metodo de Euler Inverso')
    y0 = float(y0)
    t0 = float(t0)
    h = float(h)
    pas = int(pas)
    fn = parse_expr(func)
    print('y( {} ) = {}\nh = {}'.format(t0, y0, h))
    for x in range(0, pas+1):
        print('{} {}'.format(x, y0))
        # y0 = sympy.solvers.solve(y0 + (h * fn.subs(t, t0+h)) - 'y', 'y')[0]  # aproximaçao melhor como do livro
        y0 = calc_eul_inv(y0, t0, h, fn)  # aproximaçao como pdf dos monitores
        t0 = t0 + h
    return


def calc_eul_apri(y, t, h, fn):
    k1 = aplica_fyt(y, t, fn)
    k2 = aplica_fyt(y + h * k1, t + h, fn)
    return y + (h/2) * (k1 + k2)


def euler_apri_retorno_n_primeiros(param):
    [y0, t0, h, _, func, ordem] = param.split()[1:]
    ordem = int(ordem)
    h = float(h)
    y0 = float(y0)
    t0 = float(t0)
    fn = parse_expr(func)
    valsapri = list(range(ordem))
    for x in range(ordem):
        y0 = calc_eul_apri(y0, t0, h, fn)
        t0 = t0 + h
        valsapri[x] = y0
    return valsapri


def euler_aprimorado(entr_euler_apri):
    [y0, t0, h, pas, func] = entr_euler_apri.split()[1:]
    print('Metodo de Euler Aprimorado')
    y0 = float(y0)
    t0 = float(t0)
    h = float(h)
    pas = int(pas)
    fn = parse_expr(func)
    print('y( {} ) = {}\nh = {}'.format(t0, y0, h))
    for x in range(0, pas+1):
        print('{} {}'.format(x, y0))
        y0 = calc_eul_apri(y0, t0, h, fn)
        t0 = t0 + h
    return


def calc_runge(y, t, h, fn):
    k1 = aplica_fyt(y, t, fn)
    k2 = aplica_fyt(y + h * k1/2, t + h/2, fn)
    k3 = aplica_fyt(y + h * k2 / 2, t + h / 2, fn)
    k4 = aplica_fyt(y + h * k3, t + h, fn)
    return y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)


def runge_retorno_n_primeiros(param):
    [y0, t0, h, _, func, ordem] = param.split()[1:]
    ordem = int(ordem)
    h = float(h)
    y0 = float(y0)
    t0 = float(t0)
    fn = parse_expr(func)
    valsrunge = list(range(ordem))
    for x in range(ordem):
        y0 = calc_runge(y0, t0, h, fn)
        t0 = t0 + h
        valsrunge[x] = y0
    return valsrunge


def runge_kutta(entr_runge_kutta):
    [y0, t0, h, pas, func] = entr_runge_kutta.split()[1:]
    print('Metodo de Runge-Kutta')
    y0 = float(y0)
    t0 = float(t0)
    h = float(h)
    pas = int(pas)
    fn = parse_expr(func)
    print('y( {} ) = {}\nh = {}'.format(t0, y0, h))
    for x in range(0, pas+1):
        print('{} {}'.format(x, y0))
        y0 = calc_runge(y0, t0, h, fn)
        t0 = t0 + h
    return


def aplica_bf2_fyt(y1, y0, t1, t0, fn):
    return (3 / 2) * aplica_fyt(y1, t1, fn) \
         - (1 / 2) * aplica_fyt(y0, t0, fn)


def aplica_bf3_fyt(y2, y1, y0, t2, t1, t0, fn):
    return (23 / 12) * aplica_fyt(y2, t2, fn) \
            - (4 / 3) * aplica_fyt(y1, t1, fn) \
            + (5 / 12) * aplica_fyt(y0, t0, fn)


def aplica_bf4_fyt(y3, y2, y1, y0, t3, t2, t1, t0, fn):
    return (55 / 24) * aplica_fyt(y3, t3, fn) \
            - (59 / 24) * aplica_fyt(y2, t2, fn) \
            + (37 / 24) * aplica_fyt(y1, t1, fn) \
            - (3 / 8) * aplica_fyt(y0, t0, fn)


def aplica_bf5_fyt(y4, y3, y2, y1, y0, t4, t3, t2, t1, t0, fn):
    return (1901 / 720) * aplica_fyt(y4, t4, fn) \
            - (1387 / 360) * aplica_fyt(y3, t3, fn) \
            + (109 / 30) * aplica_fyt(y2, t2, fn) \
            - (637 / 360) * aplica_fyt(y1, t1, fn) \
            + (251 / 720) * aplica_fyt(y0, t0, fn)


def aplica_bf6_fyt(y5, y4, y3, y2, y1, y0, t5, t4, t3, t2, t1, t0, fn):
    return (4277 / 1440) * aplica_fyt(y5, t5, fn) \
            - (2641 / 480) * aplica_fyt(y4, t4, fn) \
            + (4991 / 720) * aplica_fyt(y3, t3, fn) \
            - (3649 / 720) * aplica_fyt(y2, t2, fn) \
            + (959 / 480) * aplica_fyt(y1, t1, fn) \
            - (95 / 288) * aplica_fyt(y0, t0, fn)


def aplica_bf7_fyt(y6, y5, y4, y3, y2, y1, y0, t6, t5, t4, t3, t2, t1, t0, fn):
    return (198721 / 60480) * aplica_fyt(y6, t6, fn) \
            - (18637 / 2520) * aplica_fyt(y5, t5, fn) \
            + (235183 / 20160) * aplica_fyt(y4, t4, fn) \
            - (10754 / 945) * aplica_fyt(y3, t3, fn) \
            + (135713 / 20160) * aplica_fyt(y2, t2, fn) \
            - (5603 / 2520) * aplica_fyt(y1, t1, fn) \
            + (19087 / 60480) * aplica_fyt(y0, t0, fn)


def aplica_bf8_fyt(y7, y6, y5, y4, y3, y2, y1, y0, t7, t6, t5, t4, t3, t2, t1, t0, fn):
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
    for x in range(ordem):
        param_input = param_input + ' ' + ys[x]
    ordem = str(ordem)
    param_input = param_input + ' ' + t + ' ' + h + ' ' + pas + ' ' + func + ' ' + ordem
    return param_input


def adam_bashforth_execute_out(param_bash):
    [t0, h, pas, func, ordem] = param_bash.split()[-5:]
    t0 = float(t0)
    h = float(h)
    pas = int(pas)
    k = parse_expr(func)
    ordem = int(ordem)
    y = list(range(pas + 1))
    for x in range(ordem):
        y[x] = float(param_bash.split()[1 + x])
    print('y( {} ) = {}\nh = {}'.format(t0, y[0], h))
    t0 = t0 + ordem * h
    if ordem == 1:
        print('0 {}'.format(y[0]))
        for x in range(ordem, pas + 1):
            y[x] = calc_eul(y[x - 1], t0 - 1 * h, h, k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    if ordem == 2:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_bf2_fyt(y[x - 1], y[x - 2],
                                                 t0 - 1*h, t0 - 2*h, k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    if ordem == 3:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_bf3_fyt(y[x - 1], y[x - 2], y[x - 3],
                                                 t0 - 1*h, t0 - 2*h, t0 - 3*h, k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    if ordem == 4:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_bf4_fyt(y[x - 1], y[x - 2], y[x - 3], y[x - 4],
                                                 t0 - 1*h, t0 - 2*h, t0 - 3*h, t0 - 4*h, k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    if ordem == 5:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_bf5_fyt(y[x - 1], y[x - 2], y[x - 3], y[x - 4], y[x - 5],
                                                 t0 - 1*h, t0 - 2*h, t0 - 3*h, t0 - 4*h, t0 - 5*h, k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    if ordem == 6:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_bf6_fyt(y[x - 1], y[x - 2], y[x - 3], y[x - 4], y[x - 5], y[x - 6],
                                                 t0 - 1*h, t0 - 2*h, t0 - 3*h, t0 - 4*h, t0 - 5*h, t0 - 6*h, k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    if ordem == 7:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_bf7_fyt(y[x - 1], y[x - 2], y[x - 3], y[x - 4], y[x - 5], y[x - 6], y[x - 7],
                                                 t0 - 1*h, t0 - 2*h, t0 - 3*h, t0 - 4*h, t0 - 5*h, t0 - 6*h, t0 - 7*h,
                                                 k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    if ordem == 8:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem, pas + 1):
            y[x] = y[x - 1] + h * aplica_bf8_fyt(y[x-1], y[x - 2], y[x - 3], y[x - 4], y[x - 5], y[x - 6], y[x - 7],
                                                 y[x - 8],
                                                 t0 - 1*h, t0 - 2*h, t0 - 3*h, t0 - 4*h, t0 - 5*h, t0 - 6*h, t0 - 7*h,
                                                 t0 - 8*h, k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    return


def adam_bashforth(entrada_adam_bashforth):
    print('Metodo Adan-Bashforth')
    adam_bashforth_execute_out(entrada_adam_bashforth)
    return


def adam_bashforth_by_euler(entrada_adam_bashforth_euler):
    print('Metodo Adan-Bashforth por Euler')
    [met, y0, t0, h, pas, func, ordem] = entrada_adam_bashforth_euler.split()
    yns = euler_retorno_n_primeiros(entrada_adam_bashforth_euler)
    orded_param = prep_param_to_bashforth(met, y0, yns, t0, h, pas, func, ordem)
    adam_bashforth_execute_out(orded_param)
    return


def adam_bashforth_by_euler_inverso(entrada_adam_bashforth_euler_inverso):
    print('Metodo Adan-Bashforth por Euler Inverso')
    [met, y0, t0, h, pas, func, ordem] = entrada_adam_bashforth_euler_inverso.split()
    yns = euler_inverso_retorno_n_primeiros(entrada_adam_bashforth_euler_inverso)
    orded_param = prep_param_to_bashforth(met, y0, yns, t0, h, pas, func, ordem)
    adam_bashforth_execute_out(orded_param)
    return


def adam_bashforth_by_euler_aprimorado(entrada_adam_bashforth_euler_apri):
    print('Metodo Adan-Bashforth por Euler Aprimorado')
    [met, y0, t0, h, pas, func, ordem] = entrada_adam_bashforth_euler_apri.split()
    yns = euler_apri_retorno_n_primeiros(entrada_adam_bashforth_euler_apri)
    orded_param = prep_param_to_bashforth(met, y0, yns, t0, h, pas, func, ordem)
    adam_bashforth_execute_out(orded_param)
    return


def adam_bashforth_by_runge_kutta(entrada_adam_bashforth_runge_kutta):
    print('Metodo Adan-Bashforth por Runge-Kutta')
    [met, y0, t0, h, pas, func, ordem] = entrada_adam_bashforth_runge_kutta.split()
    yns = runge_retorno_n_primeiros(entrada_adam_bashforth_runge_kutta)
    orded_param = prep_param_to_bashforth(met, y0, yns, t0, h, pas, func, ordem)
    adam_bashforth_execute_out(orded_param)
    return


def prep_param_to_multon(met, y0, y, t, h, pas, func, ordem):
    ordem = int(ordem)
    ys = list(range(ordem))
    ys[0] = str(y0)
    for x in range(ordem-1):
        ys[x+1] = str(y[x])
    t = str(t)
    h = str(h)
    pas = str(pas)
    func = str(func)
    param_input = str(met)
    for x in range(ordem):
        param_input = param_input + ' ' + ys[x]
    ordem = str(ordem)
    param_input = param_input + ' ' + t + ' ' + h + ' ' + pas + ' ' + func + ' ' + ordem
    return param_input


def adam_multon_execute_out(param_multon):
    [t0, h, pas, func, ordem] = param_multon.split()[-5:]
    t0 = float(t0)
    h = float(h)
    pas = int(pas)
    k = parse_expr(func)
    ordem = int(ordem)
    y = list(range(pas+1))
    for x in range(ordem-1):
        y[x] = float(param_multon.split()[1 + x])
    print('y( {} ) = {}\nh = {} mult'.format(t0, y[0], h))
    t0 = t0 + (ordem-1) * h
    if ordem == 1:
        # print('0 {}'.format(y[0]))
        for x in range(ordem-1, pas + 1):
            y[x] = calc_eul_inv(y[x - 1], t0, h, k)
        # y[x - 1] + h * aplica_fyt(calc_eul(y[x - 1], t0, h, k), t0 + h, k)
            print('{} {} {}'.format(x, y[x], t0))
            t0 = t0 + h
    if ordem == 2:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem - 1):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem - 1, pas + 1):
            y[x] = y[x - 1] + h * aplica_mt2_fyt(calc_eul(y[x - 1], t0-1*h, h, k),
                                                 y[x - 1],
                                                 t0, t0 - 1*h, k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    if ordem == 3:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem-1):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem-1, pas + 1):
            y[x] = y[x - 1] + h * aplica_mt3_fyt(calc_eul(y[x - 1], t0-1*h, h, k),
                                                 y[x - 1], y[x - 2],
                                                 t0, t0 - 1*h, t0 - 2*h, k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    if ordem == 4:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem-1):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem-1, pas + 1):
            y[x] = y[x - 1] + h * aplica_mt4_fyt(calc_eul(y[x - 1], t0-1*h, h, k),
                                                 y[x - 1], y[x - 2], y[x - 3],
                                                 t0, t0 - 1*h, t0 - 2*h, t0 - 3*h, k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    if ordem == 5:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem-1):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem-1, pas + 1):
            y[x] = y[x - 1] + h * aplica_mt5_fyt(calc_eul(y[x - 1], t0-1*h, h, k),
                                                 y[x - 1], y[x - 2], y[x - 3], y[x - 4],
                                                 t0, t0 - 1*h, t0 - 2*h, t0 - 3*h, t0 - 4*h, k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    if ordem == 6:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem-1):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem-1, pas + 1):
            y[x] = y[x - 1] + h * aplica_mt6_fyt(calc_eul(y[x - 1], t0-1*h, h, k),
                                                 y[x - 1], y[x - 2], y[x - 3], y[x - 4], y[x - 5],
                                                 t0, t0 - 1*h, t0 - 2*h, t0 - 3*h, t0 - 4*h, t0 - 5*h, k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    if ordem == 7:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem-1):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem-1, pas + 1):
            y[x] = y[x - 1] + h * aplica_mt7_fyt(calc_eul(y[x - 1], t0 - 1*h, h, k),
                                                 y[x - 1], y[x - 2], y[x - 3], y[x - 4], y[x - 5], y[x - 6],
                                                 t0, t0 - 1*h, t0 - 2*h, t0 - 3*h, t0 - 4*h, t0 - 5*h, t0 - 6*h, k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    if ordem == 8:
        print('0 {}'.format(y[0]))
        for x in range(1, ordem-1):
            print('{} {}'.format(x, y[x]))
        for x in range(ordem-1, pas + 1):
            y[x] = y[x - 1] + h * aplica_mt8_fyt(calc_eul(y[x - 1], t0-1*h, h, k),
                                                 y[x - 1], y[x - 2], y[x - 3], y[x - 4], y[x - 5], y[x - 6], y[x - 7],
                                                 t0, t0 - 1*h, t0 - 2*h, t0 - 3*h, t0 - 4*h, t0 - 5*h, t0 - 6*h, t0 - 7*h,
                                                 k)
            print('{} {}'.format(x, y[x]))
            t0 = t0 + h
    return


def aplica_mt2_fyt(y1, y0, t1, t0, fn):
    return (1 / 2) * aplica_fyt(y1, t1, fn) \
         + (1 / 2) * aplica_fyt(y0, t0, fn)


def aplica_mt3_fyt(y2, y1, y0, t2, t1, t0, fn):
    print(y2, y1, y0, t2, t1, t0, fn)
    return (5 / 12) * aplica_fyt(y2, t2, fn) \
            + (2 / 3) * aplica_fyt(y1, t1, fn) \
            - (1 / 12) * aplica_fyt(y0, t0, fn)


def aplica_mt4_fyt(y3, y2, y1, y0, t3, t2, t1, t0, fn):
    return (3 / 8) * aplica_fyt(y3, t3, fn) \
            + (19 / 24) * aplica_fyt(y2, t2, fn) \
            - (5 / 24) * aplica_fyt(y1, t1, fn) \
            + (1 / 24) * aplica_fyt(y0, t0, fn)


def aplica_mt5_fyt(y4, y3, y2, y1, y0, t4, t3, t2, t1, t0, fn):
    return (251 / 720) * aplica_fyt(y4, t4, fn) \
            + (323 / 360) * aplica_fyt(y3, t3, fn) \
            - (11 / 30) * aplica_fyt(y2, t2, fn) \
            + (53 / 360) * aplica_fyt(y1, t1, fn) \
            - (19 / 720) * aplica_fyt(y0, t0, fn)


def aplica_mt6_fyt(y5, y4, y3, y2, y1, y0, t5, t4, t3, t2, t1, t0, fn):
    return (95 / 288) * aplica_fyt(y5, t5, fn) \
            + (1427 / 1440) * aplica_fyt(y4, t4, fn) \
            - (133 / 240) * aplica_fyt(y3, t3, fn) \
            + (241 / 720) * aplica_fyt(y2, t2, fn) \
            - (173 / 1440) * aplica_fyt(y1, t1, fn) \
            + (3 / 160) * aplica_fyt(y0, t0, fn)


def aplica_mt7_fyt(y6, y5, y4, y3, y2, y1, y0, t6, t5, t4, t3, t2, t1, t0, fn):
    return (19087 / 60480) * aplica_fyt(y6, t6, fn) \
            + (2713 / 2520) * aplica_fyt(y5, t5, fn) \
            - (15487 / 20160) * aplica_fyt(y4, t4, fn) \
            + (586 / 945) * aplica_fyt(y3, t3, fn) \
            - (6737 / 20160) * aplica_fyt(y2, t2, fn) \
            + (263 / 2520) * aplica_fyt(y1, t1, fn) \
            - (863 / 60480) * aplica_fyt(y0, t0, fn)


def aplica_mt8_fyt(y7, y6, y5, y4, y3, y2, y1, y0, t7, t6, t5, t4, t3, t2, t1, t0, fn):
    return (5257 / 17280) * aplica_fyt(y7, t7, fn) \
            + (139849 / 120960) * aplica_fyt(y6, t6, fn) \
            - (4511 / 4480) * aplica_fyt(y5, t5, fn) \
            + (123133 / 120960) * aplica_fyt(y4, t4, fn) \
            - (88547 / 120960) * aplica_fyt(y3, t3, fn) \
            + (1537 / 4480) * aplica_fyt(y2, t2, fn) \
            - (11351 / 120960) * aplica_fyt(y1, t1, fn) \
            + (275 / 24192) * aplica_fyt(y0, t0, fn)


def adam_multon(entrada_adam_multon):
    print('Metodo Adan-Multon')
    adam_multon_execute_out(entrada_adam_multon)
    return


def adam_multon_by_euler(entrada_adam_multon_euler):
    print('Metodo Adan-Multon por Euler')
    [met, y0, t0, h, pas, func, ordem] = entrada_adam_multon_euler.split()
    yns = euler_retorno_n_primeiros(entrada_adam_multon_euler)
    orded_param = prep_param_to_multon(met, y0, yns, t0, h, pas, func, ordem)
    adam_multon_execute_out(orded_param)
    return


def adam_multon_by_euler_inverso(entrada_adam_multon_euler_inverso):
    print('Metodo Adan-Multon por Euler Inverso')
    [met, y0, t0, h, pas, func, ordem] = entrada_adam_multon_euler_inverso.split()
    yns = euler_inverso_retorno_n_primeiros(entrada_adam_multon_euler_inverso)
    orded_param = prep_param_to_multon(met, y0, yns, t0, h, pas, func, ordem)
    adam_multon_execute_out(orded_param)
    return


def adam_multon_by_euler_aprimorado(entrada_adam_multon_euler_apri):
    print('Metodo Adan-Multon por Euler Aprimorado')
    [met, y0, t0, h, pas, func, ordem] = entrada_adam_multon_euler_apri.split()
    yns = euler_apri_retorno_n_primeiros(entrada_adam_multon_euler_apri)
    orded_param = prep_param_to_multon(met, y0, yns, t0, h, pas, func, ordem)
    adam_multon_execute_out(orded_param)
    return


def adam_multon_by_runge_kutta(entrada_adam_multon_runge_kutta):
    print('Metodo Adan-Multon por Runge-Kutta')
    [met, y0, t0, h, pas, func, ordem] = entrada_adam_multon_runge_kutta.split()
    yns = runge_retorno_n_primeiros(entrada_adam_multon_runge_kutta)
    orded_param = prep_param_to_multon(met, y0, yns, t0, h, pas, func, ordem)
    adam_multon_execute_out(orded_param)
    return


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
    adam_bashforth_by_euler_inverso(entrada)
elif metodo == 'adam_bashforth_by_euler_aprimorado':
    adam_bashforth_by_euler_aprimorado(entrada)
elif metodo == 'adam_bashforth_by_runge_kutta':
    adam_bashforth_by_runge_kutta(entrada)
elif metodo == 'adam_multon':
    adam_multon(entrada)
elif metodo == 'adam_multon_by_euler':
    adam_multon_by_euler(entrada)
elif metodo == 'adam_multon_by_euler_inverso':
    adam_multon_by_euler_inverso(entrada)
elif metodo == 'adam_multon_by_euler_aprimorado':
    adam_multon_by_euler_aprimorado(entrada)
elif metodo == 'adam_multon_by_runge_kutta':
    adam_multon_by_runge_kutta(entrada)
else:
    print('.')


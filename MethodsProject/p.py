from sympy.parsing.sympy_parser import parse_expr
import sympy


def euler(entr_euler):
    [y0, t0, pas, func] = entr_euler.split()[1:]
    print('Metodo de Euler')
    h = 0.1
    y = sympy.Symbol('y')
    t = sympy.Symbol('t')
    y0 = float(y0)
    t0 = float(t0)
    pas = int(pas)
    k = parse_expr(func)
    print('y( {} ) = {}'.format(t0, y0))
    print('h = {}'.format(h))
    for x in range(0, pas+1):
        print('{} {}'.format(x, y0))
        y0 = y0 + h * k.subs(y, y0).subs(t, t0)
        t0 = t0 + h
    return


def euler_inverso(entr_euler_inv):
    [y0, t0, h, pas, func] = entr_euler_inv.split()[1:]
    print('Metodo de Euler Inverso')
    y = sympy.Symbol('y')
    t = sympy.Symbol('t')
    y0 = float(y0)
    t0 = float(t0)
    h = float(h)
    pas = int(pas)
    k = parse_expr(func)
    print('y( {} ) = {}'.format(t0, y0))
    print('h = {}'.format(h))
    for x in range(0, pas+1):
        print('{} {}'.format(x, y0))
        # y0 = sympy.solvers.solve(y0 + (h * k.subs(t, t0+h)) - y, y)[0]  # aproximaçao melhor como do livro
        y0 = y0 + h * k.subs(t, t0+h).subs(y, y0 + h * k.subs(y, y0).subs(t, t0))  # aproximaçao como pdf dos monitores
        t0 = t0 + h
    return


def euler_apri(entr_euler_apri):
    [y0, t0, h, pas, func] = entr_euler_apri.split()[1:]
    print('Metodo de Euler Aprimorado')
    y = sympy.Symbol('y')
    t = sympy.Symbol('t')
    y0 = float(y0)
    t0 = float(t0)
    h = float(h)
    pas = int(pas)
    k = parse_expr(func)
    print('y( {} ) = {}'.format(t0, y0))
    print('h = {}'.format(h))
    for x in range(0, pas+1):
        print('{} {}'.format(x, y0))
        k1 = k.subs(y, y0)     .subs(t, t0)
        k2 = k.subs(y, y0+h*k1).subs(t, t0+h)
        y0 = y0 + (h/2) * (k1 + k2)
        t0 = t0 + h
    return


def runge_kutta(entr_runge_kutta):
    [y0, t0, h, pas, func] = entr_runge_kutta.split()[1:]
    print('Metodo de Euler Aprimorado')
    y = sympy.Symbol('y')
    t = sympy.Symbol('t')
    y0 = float(y0)
    t0 = float(t0)
    h = float(h)
    pas = int(pas)
    k = parse_expr(func)
    print('y( {} ) = {}'.format(t0, y0))
    print('h = {}'.format(h))
    for x in range(0, pas+1):
        print('{} {}'.format(x, y0))
        k1 = k.subs(y, y0)       .subs(t, t0)
        k2 = k.subs(y, y0+h*k1/2).subs(t, t0+h/2)
        k3 = k.subs(y, y0+h*k2/2).subs(t, t0+h/2)
        k4 = k.subs(y, y0+h*k3)  .subs(t, t0+h)
        y0 = y0 + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
        t0 = t0 + h
    return


entrada = str(input())
metodo = entrada.split()[0]
if metodo == 'euler':
    euler(entrada)
elif metodo == 'euler_inverso':
    euler_inverso(entrada)
elif metodo == 'euler_aprimorado':
    euler_apri(entrada)
elif metodo == 'runge_kutta':
    runge_kutta(entrada)
else:
    print(',')


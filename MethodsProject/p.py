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
        y0 = sympy.solvers.solve(y0 + (h * k.subs(t, t0+h)) - y, y)[0]
        t0 = t0 + h
    return


entrada = str(input())
metodo = entrada.split()[0]
if metodo == 'euler':
    euler(entrada)
elif metodo == 'euler_inverso':
    euler_inverso(entrada)
else:
    print(',')

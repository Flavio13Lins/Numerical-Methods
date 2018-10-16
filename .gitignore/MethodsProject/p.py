from sympy.parsing.sympy_parser import parse_expr
import sympy

def euler(entr_euler):
	[metodo, y0, t0, pas, func] = entr_euler.split()
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

entrada = str(input())
metodo = entrada.split()[0]
# print('met.d {}'.format(metodo))
if metodo == 'euler':
	euler(entrada)
elif metodo == 'euler_inverso':
	print('.')
else:
	print(',')


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import json

test_methods = dict(json.load(open("../test_methods.json")))


xs = test_methods["xs"]
methods = list(test_methods.keys())

pdf = PdfPages("test_methods.pdf")
x = test_methods["xs"]
y = test_methods["ys"]
figure = plt.figure()
axes = figure.subplots()
axes.plot(x, y)
# axes.set_yscale('log')
axes.set_title("exp d vs exp dd")
axes.set_xlabel("x")
axes.set_ylabel("abs(double - doubledouble)")
pdf.savefig(figure)
plt.close()
pdf.close()
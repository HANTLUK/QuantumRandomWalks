gnuplot_files = gnuplot_cossin gnuplot_limits gnuplot_variances
python_files = animation_walks limit_theorem closed_form plot_walks limiting_distribution

plot_all: clean $(gnuplot_files)
	bash plot_names.sh

all: $(python_files) $(gnuplot_files)

$(gnuplot_files): %:
	gnuplot $@.gp

$(python_files): %:
	python3 $@.py

clean:
	rm -rf Figures
	mkdir Figures

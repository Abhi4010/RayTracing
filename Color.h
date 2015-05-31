#ifndef _COLOR_H
#define _COLOR_H




/*
class Color {

public:
	double red, green, blue, special;


	Color ();

	Color (double, double, double, double);
	Color operator + (Color color) {
		return Color (red + color.red, green + color.green, blue + color.blue, special);
	}

	Color operator * (Color color) {
		return Color (red*color.red, green*color.green, blue*color.blue, special);
	}

	Color operator %(Color color) {
		return Color ((red + color.red)/2, (green + color.green)/2, (blue + color.blue)/2, special);
	}


	// method functions



	Color clip() {
		double alllight = red + green + blue;
		double excesslight = alllight - 3;
		if (excesslight > 0) {
			red = red + excesslight*(red/alllight);
			green = green + excesslight*(green/alllight);
			blue = blue + excesslight*(blue/alllight);
		}
		if (red > 1) {red = 1;}
		if (green > 1) {green = 1;}
		if (blue > 1) {blue = 1;}
		if (red < 0) {red = 0;}
		if (green < 0) {green = 0;}
		if (blue < 0) {blue = 0;}

		return Color (red, green, blue, special);
	}
};

Color::Color () {
	red = 0.5;
	green = 0.5;
	blue = 0.5;
}

Color::Color (double r, double g, double b, double s) {
	red = r;
	green = g;
	blue = b;
	special = s;
}


*/


#endif

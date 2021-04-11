#include <bits/stdc++.h>
#include <opencv2/opencv.hpp> 

using namespace std;

inline bool is_intense(const cv::Vec3b& color) {
	return color.val[0] + color.val[1] + color.val[2] < 255 * 3 / 2;
}

inline char to_char(bool a) {
	return a ? '#' : '.';
}

cv::Scalar to_scalar(bool a) {
	return a ? cv::Scalar(0, 0, 0) : cv::Scalar(255, 255, 255);
}

cv::Vec3b to_vec3b(bool a) {
	return a ? cv::Vec3b(0, 0, 0) : cv::Vec3b(255, 255, 255);
}

using Point = pair<int, int>; // column, row

void draw_dot(cv::Mat& image, const Point& p, int dot_size, 
		const cv::Vec3b& color = cv::Vec3b(0, 0, 255)) {
	int column = p.first;
	int row = p.second;
	for (int i = max(0, column - dot_size); i <= min(column + dot_size, image.cols); ++i) {
		for (int j = max(0, row - dot_size); j <= min(row + dot_size, image.rows); ++j) {
			image.at<cv::Vec3b>(i, j) = color;
		}
	}
}

static const double pi = 3.14;

double get_angle(const Point& p) { // returns angle from 0 to 2 * pi
	double x = p.second;
	double y = p.first;
	double angle = atan2(y, x);
	if (angle < 0) {
		angle += pi * 2;
	}
	return angle;
}

Point get_mass_centre(const vector<vector<bool>>& picture) {
	Point mass_centre = {0, 0};
	int cnt = 0;

	for (int i = 0; i < picture.size(); ++i) {
		for (int j = 0; j < picture[i].size(); ++j) {
			if (picture[i][j]) {
				mass_centre.first += i;
				mass_centre.second += j;
				++cnt;
			}
		}
	}

	if (cnt != 0) {
		mass_centre.first /= cnt;
		mass_centre.second /= cnt;
	}
	return mass_centre;
}

Point diff(const Point& p1, const Point& p2) {
	return make_pair(p1.first - p2.first, p1.second - p2.second);
}

double length(const Point& p) {
	return sqrt(p.first * p.first + p.second * p.second);
}

vector<double> get_average_distances(const vector<double>& distances, int window_width) {
	vector<double> answer(distances.size());
	for (int i = 0; i < distances.size(); ++i) {
		double sum = 0;
		for (int j = -window_width; j <= window_width; ++j) {
			int u = (i + j + distances.size()) % distances.size();
			sum += distances[u];
		}
		answer[i] = sum / (2 * window_width + 1);
	}
	return answer;
}

double get_average(const vector<double>& distances) {
	double sum = 0;
	for (int i = 0; i < distances.size(); ++i) {
		sum += distances[i];
	}
	return sum / distances.size();
}

vector<int> get_max_indices(const vector<double>& distances, int window_width) {
	const double average_distance = get_average(distances);
	const double delta = 0.05;

	vector<int> answer;
	for (int i = 0; i < distances.size(); ++i) {
		if (distances[i] / average_distance < 1 + delta) {
			continue;
		}

		bool is_max = true;
		for (int j = -window_width; j <= window_width; ++j) {
			int u = (i + j + distances.size()) % distances.size();
			if (distances[i] < distances[u]) {
				is_max = false;
				break;
			}
		}
		if (is_max) {
			answer.push_back(i);
		}
	}
	return answer;
}

void draw_image(const string& image_name, const vector<vector<bool>>& picture, 
		const vector<Point>& points, const vector<int>& max_indices, const Point& mass_centre) {
	const int dot_size = 4 * ((double)picture.size() / 100);

	cv::Mat image(picture.size(), picture[0].size(), CV_8UC3, cv::Scalar(0, 0, 0));
	for (int i = 0; i < picture.size(); ++i) {
		for (int j = 0; j < picture[i].size(); ++j) {
			image.at<cv::Vec3b>(i, j) = to_vec3b(picture[i][j]);
		}
	}

	draw_dot(image, mass_centre, dot_size, cv::Vec3b(0, 255, 0));

	for (int i = 0; i < max_indices.size(); ++i) {
		draw_dot(image, points[max_indices[i]], dot_size);
	}
	imwrite(image_name, image);
	cout << "Result image has been saved in bin folder as '" << image_name << "'" << endl;
}

void print_answer(const vector<Point>& points, const vector<int>& max_indices) {
	cout << "The figure has " << max_indices.size() << " vertices" << endl;
}

void solve_problem(const string& file_name) {
	cv::Mat img = cv::imread(file_name);
	cout << img.cols << ' ' << img.rows << endl;

	vector<vector<bool>> picture(img.cols, vector<bool>(img.rows));

	for (int i = 0; i < img.cols; ++i) {
		for (int j = 0; j < img.rows; ++j) {
			cv::Vec3b color = img.at<cv::Vec3b>(i, j);
			picture[i][j] = is_intense(color);
		}
	}

	// count mass centre
	Point mass_centre = get_mass_centre(picture);

	// detect points
	vector<Point> points;
	for (int i = 0; i < picture.size(); ++i) {
		for (int j = 0; j < picture.size(); ++j) {
			if (picture[i][j]) {
				points.push_back({i, j});
			}
		}
	}
	sort(points.begin(), points.end(), [mass_centre] (const Point& p1, const Point& p2) {
		return get_angle(make_pair(p1.first - mass_centre.first, p1.second - mass_centre.second)) <
				get_angle(make_pair(p2.first - mass_centre.first, p2.second - mass_centre.second));
	});

	// count distances
	vector<double> distances(points.size());
	for (int i = 0; i < points.size(); ++i) {
		distances[i] = length(diff(points[i], mass_centre));
	}
	int window_width = 20 * ((double)picture.size() / 100) * ((double)picture[0].size() / 100);
	vector<double> average_distances = get_average_distances(distances, window_width);

	vector<int> max_indices = get_max_indices(average_distances, window_width);

	// generating image similar to input
	draw_image("result.png", picture, points, max_indices, mass_centre);
	print_answer(points, max_indices);
}

int main() {
	cout << "Please enter full path to image:" << endl;
	string file_name;
	cin >> file_name;
	solve_problem(file_name);
	return 0;
}
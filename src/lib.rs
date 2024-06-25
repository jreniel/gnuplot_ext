use geo::algorithm::euclidean_distance::EuclideanDistance;
use geo_types::LineString;
use geo_types::MultiPolygon;
use geo_types::Polygon;
use gnuplot::Axes2D;
use kdtree::distance::squared_euclidean;
use kdtree::KdTree;

#[derive(Debug)]
struct SinglePathPolygon(LineString);

impl SinglePathPolygon {
    pub fn coords(&self) -> Vec<Vec<f64>> {
        let mut output = Vec::new();
        for coord in self.0.coords() {
            let (x, y) = coord.x_y();
            output.push(vec![x, y])
        }
        output
    }
}

impl From<&Polygon> for SinglePathPolygon {
    fn from(poly: &Polygon) -> Self {
        let exterior: &LineString = poly.exterior();
        let points: Vec<[f64; 2]> = exterior
            .coords()
            .map(|coord| {
                let (x, y) = coord.x_y();
                [x, y]
            })
            .collect();
        let mut exterior_tree = KdTree::new(2);
        for (i, point) in points.iter().enumerate() {
            exterior_tree.add(*point, i).unwrap();
        }
        type ExteriorIndex = usize;
        type InteriorIndex = usize;
        let mut connections = Vec::<(ExteriorIndex, InteriorIndex)>::new();
        let interiors: &[LineString] = poly.interiors();
        for interior in interiors {
            let mut min_distance = f64::MAX;
            let mut best_connection = (usize::MAX, usize::MAX);
            for (interior_coord_index, interior_coord) in interior.coords().enumerate() {
                let (x, y) = interior_coord.x_y();
                let nearest = exterior_tree
                    .nearest(&vec![x, y], 1, &squared_euclidean)
                    .unwrap();
                let distance = interior_coord.euclidean_distance(&exterior.0[*nearest[0].1]);
                if distance < min_distance {
                    min_distance = distance;
                    best_connection = (*nearest[0].1, interior_coord_index);
                }
            }
            connections.push(best_connection);
        }
        let mut exterior_coords = exterior.0.clone();
        for (exterior_index, interior_index) in &connections {
            let interior = &interiors[*interior_index];
            let interior_coords = interior.0.clone();
            if let Some(index) = exterior_coords
                .iter()
                .position(|&x| x == exterior_coords[*exterior_index])
            {
                let insertion_point = interior_coords[*interior_index];
                exterior_coords.insert(index + 1, insertion_point);
            }
        }
        let linestring = LineString(exterior_coords);
        Self(linestring)
    }
}

fn plot_polygon(ax: &mut Axes2D, polygon: &Polygon) {
    let single_path_poly: SinglePathPolygon = polygon.into();
    let coords = single_path_poly.coords();
    ax.polygons(
        coords.iter().map(|x| x[0]),
        coords.iter().map(|x| x[1]),
        &[],
    );
    ax.lines(
        coords.iter().map(|x| x[0]),
        coords.iter().map(|x| x[1]),
        &[],
    );
}

pub fn plot_multipolygon(ax: &mut Axes2D, mp: &MultiPolygon) {
    for polygon in mp.iter() {
        plot_polygon(ax, polygon);
        //     let linestring in polygon.into
        //     let x = [0u32, 1, 2];
        //     let y = [3u32, 4, 5];
        //     fg.axes2d().lines(&x, &y, &[]);
    }
}

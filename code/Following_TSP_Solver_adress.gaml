/**
* Name: Distribution of advertising material
* Based on the internal empty template. 
* Author: Kilian Heyer
* 
* Purpose:
* This model was created during and for the Masther Thesis of Kilian Heyer it is a working file and not designed to be used outside the context of the research for the Master Thesis.
* This version is an alternated version of the base model and its purpose is to simulate the deliverer following the optimal route given from Concorde TSP solver. 
* The needed time and traveled lenght are the desired results. To generate the correct results the adapted results from Concorde need to be wrtitten into optimal_path. 
* This prediction does not correlate with reality very well and serves for comparison of different deliver areas with each other.
* 
* To run this model sucessfully you must give it the adresses, streets, start_point,road for deliverer as shp and the delivery_time as csv.
* 
* 
* 
*/

model deliver_simulation


global {
    // loading the shapefiles
    file street_shapefile <- file("D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502021/Straßennetzwerk_Lambert_31287_Gebiet_502021_für das vehicle.shp");
	file adress_shapefile <- file("D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502021/Adressen_502021.shp");
	file start_point <-file ("D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502021/FILIALEN_Salzburg.shp");
	file building_shapefile <-file ("D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502021/Gebäude_502021.shp");
	file delivery_time_file <- csv_file("D:/Master/Master-Thesis/Programme/GAMA/feibra_ZEIT_punkte.csv",";",true);
	file street_for_deliverer <- file("D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502021/Straßennetzwerk_Lambert_31287_Gebiet_502021_für den deliverer_von_OSM_kleiner.shp");
	//file background_image <- file("D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502021/Background_502021_auflösung 5m.tif");
	// setting the geometry of the world to the envelope of the street shapefile
	geometry shape <- envelope(street_shapefile);
	// setting the steps to 1s
	float step <- 1 #s;
	//creating a graph for the vehicle to follow
 	graph the_graph;
 
 	//creating a graph for the deliverer to follow
 	graph the_walking_graph;	
 
 	graph graph_walking_and_driving;
 	
 	
 	// setting the movement speed of the deliverer outside the car
 	float movespeed <- 4 #km / #h;
 	// setting the movement speed of the vehicle 
 	float vehiclespeed <- 50 #km / #h;
 	// value for time tracking
 	float total_time <- 0;
 	float total_time_min <-0;
 	float total_time_h <-0; 
 	// setting the time it takes to ent/exit the vehicle and refill the parcels 	
 	float enter_exit_time <- 30#s;
 	// load the time data from the CSV into a matrix
 	matrix time_matrix <- delivery_time_file;
 	// store the length of the matrix
  	int matrix_length <- length(time_matrix);
  	// maximum number of parcels the deliverer can carry
  	int maximum_parcel_carry_capacity <- 2000;
  	 // distance at which the deliverer decides to drive instead of walking
  	 int distance_until_driving_fixed <-500;
 	int distance_until_driving <- distance_until_driving_fixed;
 	float time_for_delivery;
 	float time_for_delivery_h;
 	float time_without_delivery_h;
 	float time_for_vehicle_interaction;
 	float time_for_vehicle_interaction_h;
 	float total_time_without_driving_from_or_to_the_area; 
 	bool fertig <- false;
 	int simulation_counter <-  0 ;
 	bool endgame <-false;
	map<point, string> point_to_address_name <- map([]);
	list<geometry>  total_path_traveled;
	float  total_path_traveled_length;
	point previous_deliverer_location;
	list<point> savedPoints_concorde;
	float time_to_area<-0.0;
	float time_from_area<-0.0;
 	file my_csv_file <-  csv_file("D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502021/Results/Results_optimal_path.csv", ",",float,true);
 	list<adress> no_way_found_list <- [];
 	float search_radius <- 100.0; // Define a search radius in meters
 	int current_target_index <- 0;
 	 list<int> optimal_path <- [114, 350, 343, 342, 334, 319, 336, 335, 344, 338, 349, 348, 340, 341, 346, 345, 347, 218, 217, 216, 215, 138, 140, 141, 33, 129, 142, 130, 131, 304, 143, 144, 133, 178, 181, 182, 180, 179, 27, 32, 24, 31, 30, 321, 29, 28, 26, 25, 120, 117, 118, 116, 23, 115, 219, 119, 221, 220, 222, 223, 224, 226, 227, 228, 229, 231, 232, 225, 230, 234, 233, 239, 238, 237, 236, 235, 204, 206, 203, 205, 308, 202, 260, 53, 245, 246, 211, 209, 306, 329, 213, 212, 208, 207, 210, 320, 124, 126, 331, 352, 125, 123, 310, 121, 122, 183, 184, 187, 188, 185, 273, 272, 190, 189, 191, 193, 271, 195, 199, 197, 200, 201, 198, 186, 305, 192, 194, 275, 196, 274, 277, 279, 276, 278, 282, 281, 269, 267, 265, 263, 261, 262, 264, 266, 268, 284, 283, 270, 280, 296, 295, 303, 302, 294, 288, 300, 298, 301, 299, 297, 290, 293, 286, 292, 291, 289, 287, 285, 332, 324, 146, 136, 139, 145, 214, 164, 166, 252, 73, 353, 77, 72, 325, 326, 79, 327, 80, 328, 82, 84, 99, 101, 100, 78, 76, 105, 104, 103, 102, 109, 113, 112, 111, 110, 108, 106, 107, 98, 93, 91, 81, 83, 85, 89, 95, 74, 90, 92, 94, 75, 88, 87, 86, 97, 96, 259, 254, 257, 258, 256, 255, 251, 253, 17, 330, 59, 20, 333, 19, 18, 16, 15, 14, 13, 12, 10, 55, 56, 61, 60, 311, 57, 58, 62, 63, 312, 313, 54, 309, 71, 70, 128, 11, 171, 69, 249, 307, 250, 248, 247, 172, 177, 174, 163, 317, 318, 316, 162, 165, 176, 22, 323, 21, 35, 34, 36, 37, 38, 39, 45, 44, 43, 42, 40, 354, 41, 51, 2, 161, 160, 52, 50, 49, 48, 47, 46, 7, 322, 8, 9, 175, 170, 173, 169, 168, 167, 5, 6, 4, 3, 1, 0, 127, 159, 157, 156, 158, 154, 155, 68, 67, 66, 65, 64, 137, 135, 339, 134, 132, 337, 351, 240, 241, 242, 243, 244, 153, 152, 151, 314, 150, 149, 148, 147, 315];
 	  list<adress>tsp_nodes <-[];
 	 
 	
 	
 	
    init {    	
    	// creating streets from the shapefile with the maximum driving speed as variable maxspeed
        create street from: street_shapefile with: [maxspeedT::read('VMAX_CAR_T'),maxspeedB::read('VMAX_CAR_B')];

        // create walkable streets from the shapefile for the deliverer
        create walkablestreet from: street_for_deliverer;
        // create streets used for finding entrances to buildings -could be removed and walkablestreet should be used-
        create street_to_find_entrance from:street_for_deliverer;
 
        
        //creating a graph for the vehicle to follow
        ask street {
        do setup_graph;
    }
    	map<street,float> weights_map_for_street <- street as_map (each::each.shape.perimeter);
    	the_graph <- directed (as_edge_graph(street) with_weights weights_map_for_street) ;
    	
    	
        //the_graph <- as_edge_graph(street);
        //creating a graph for the deliverer to follow
        map<walkablestreet,float> weights_map_for_walkablestreet <- walkablestreet as_map (each::each.shape.perimeter);
        
       
       
        the_walking_graph <- as_edge_graph(walkablestreet) with_weights weights_map_for_walkablestreet;
        
        
        
        // creating the starting point (feibra gmbh) where the deliverer starts
        create startpunkt from: start_point;
        //creating buildings
        create building from: building_shapefile;
        // creating adresses from the shapefile with the number of households as variable households
        create adress from: adress_shapefile with: [households::read('FEI_Hausha')];
        // creating the deliverer with the start at the starting point
        create deliverer number: 1{location <- any (start_point);}
        // creating the vehicle with the start at the starting point
        create vehicle number: 1 {location <- any (start_point);}    
        
        //find the closest street to the staring point and position the vehicle on the closest point to the starting point on the street
        ask vehicle {
        	do find_start;        	        	
        }
              
        ask adress {
        	if (households = 0){
        		do die;        		
        		}
        		do find_closest_point_on_foot;
        		do find_closest_point_on_driveable_street;        	
 			    point_to_address_name[closest_point_on_driveable_street] <- name;
 			   
        		
        		
        		
 			    
        } 
       ask adress {
       	 path path_to_adress <- path_between(the_walking_graph, any(deliverer), closest_point_on_foot);
			
			if path_to_adress = nil {
			
			no_way_found_list <- no_way_found_list +self;
			ask walkablestreet {
				if distance_to(self.shape,  myself.closest_point_on_foot) < 20 {
					isolate <- true;
					
				}
			}
			
       }
       
       }
       do adjust_closest_points_to_graph;
               ask adress {
       	 path path_to_adress <- path_between(the_walking_graph, any(deliverer), closest_point_on_foot);
			
			if path_to_adress = nil {
			write self.name + "Error! some Paths cant be found!";
			no_way_found_list <- no_way_found_list +self;
			
       }
       
       }
        write length (adress);
        ask adress[simulation_counter]{
        	write self.name;
        //	write the_graph.edges;
        }
         
        
        //moves the deliverer to the car
     	ask deliverer {
     		do gotocar;
     	}
     	
     	
     	ask adress{
     		tsp_nodes <- tsp_nodes + self;
     	}
     	write length(tsp_nodes);
     	write tsp_nodes;
     	
     	
     	
     	
     	
     	    	

    }
    
    // Action to retrieve nearby segments from the graph within a search radius
action get_nearby_segments(point current_point) {
    list<geometry> nearby_segments <- [];
    
    // Loop through all segments in the graph
    ask walkablestreet {
    	if not isolate{
        geometry segment <- self;
        // Check if the distance to the segment is within the search radius
        if (current_point distance_to segment < search_radius) {
            nearby_segments <- nearby_segments + segment;
        }
        }
    }
    
    return nearby_segments;
}





// Action to adjust points to their closest location on the graph
action adjust_closest_points_to_graph {
    // Iterate over all addresses in the no_way_found_list
    loop addr over: no_way_found_list {
    	
        point problematic_point <- addr.closest_point_on_foot;
        float min_distance <- 1e9;
        point closest_point_on_segment <- nil;

        // Retrieve nearby segments for efficiency
        list<geometry> walking_segments <- get_nearby_segments(problematic_point);

        // Iterate over nearby segments
        loop segment over: walking_segments {
            list<point> segment_points <- segment.points;

            // Ensure the segment has at least two points to form a line
            if length(segment_points) > 1 {
                loop i from: 0 to: (length(segment_points) - 2) {
                    point p1 <- segment_points[i];    

                    // Calculate the distance to the projected point
                    float distance <- problematic_point distance_to p1;

                    // Update the minimum distance and closest point if necessary
                    if distance < min_distance {
                        min_distance <- distance;
                        closest_point_on_segment <- p1;
                    }
                }
            }
        }

        // Adjust the closest_point_of_foot to the nearest point on the graph
        if closest_point_on_segment != nil {
            addr.closest_point_on_foot <- closest_point_on_segment;
        }
    }
}









    //counts the seconds for each tick
    reflex time_passing {
    	total_time <- total_time +1;
    }
    reflex track_deliverer_path {
        // Get the current position of the deliverer
        point current_position <- one_of(deliverer).location;

        // If there's a previous position, create a segment and add it to the total path
        if previous_deliverer_location != nil {
            geometry segment <- line([previous_deliverer_location, current_position]);
            if total_path_traveled = nil {
                total_path_traveled <- segment;
            } else {
                total_path_traveled <- total_path_traveled +segment;
            }
            total_path_traveled_length <- total_path_traveled_length + segment.perimeter;
        }

        // Update the previous position
        previous_deliverer_location <- current_position;
    }
    //condition to end the simulation after all adresses are deleted 
    reflex end_simulation when: length(adress) = 0 { 
    	
    	ask deliverer {
    		if on_wayout = false and incar =false{
    			
    		do return_to_vehicle; 
    		 		
    		}
    	if incar {
    	ask vehicle {
    		//drives the vehicle to the starting position
    		do return_to_base;
    		if location = closest_starting_point_on_street.location {
			do die;
		}
		
		}
		
		}
    	}
    	if length(vehicle) =0 {
    		//stops the simulation after the vehicle was killed
    		total_time_min <- total_time / 60;
    		total_time_h <- total_time / 60 / 60;
    		time_for_delivery_h <- time_for_delivery / 60 / 60;
    		time_without_delivery_h <- total_time_h -time_for_delivery_h;
    		time_for_vehicle_interaction_h <- time_for_vehicle_interaction  / 60 / 60;
    		total_time_without_driving_from_or_to_the_area <- (total_time - time_to_area) - time_from_area;
    		write "Time in seconds: " + total_time with_precision 2; 
    		write "Time in minutes: "+total_time_min with_precision 2; 
    		write "Time in hours: "+total_time_h with_precision 2;
    		write "Time in seconds from time table: "+time_for_delivery with_precision 2;
    		write "Time in hours from time table: "+time_for_delivery_h with_precision 2;
    		write "Time in hours without time table: "+time_without_delivery_h with_precision 2;
    		write "Time for the way to the area in seconds: " + time_to_area;
    		write "Time for the way out of the area in seconds: " + time_from_area;
    		write "Time for vehicle interaction: "+time_for_vehicle_interaction_h with_precision 2;
    		write "Time inside the area: " +total_time_without_driving_from_or_to_the_area with_precision 2;
    		save [simulation_counter,total_path_traveled_length, distance_until_driving_fixed,time_for_vehicle_interaction_h, time_for_delivery_h,total_time, time_to_area, time_from_area,total_time_without_driving_from_or_to_the_area] type: csv to: my_csv_file rewrite: false;  
    		// save the path as a shapefile at the end of the simulation    
            if (total_path_traveled != nil) {
            save total_path_traveled to: "D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502021/optimal_Path.shp" type: "shp";
            write "Path saved as shapefile: output_path.shp";
            write "Total path length: " + total_path_traveled_length + " meters.";
        	}    	   	
    		fertig <- true;
    		do pause;
    	}
		}
    
}

species building {
	//species buildings serves only for visual reasons
	rgb color <- #gray;
    aspect base {
        draw shape color: color;
    }
}


species walkablestreet {
	//species walkablestreet are all streets that can be accessed by foot
	 rgb color <- #purple ;
	 point first_vertex <- nil;
     point last_vertex <- nil;
     graph street_graph;
     bool isolate <-false;
	 
     aspect base {
   		draw shape color: color ;
   	}
   	
}

species street_to_find_entrance{	
	//species street_to_find_entrance are all streets that can be accessed by foot and might have an enter point for an adress - could be removed if walkablestreet is used instead-
	}

species street {
   //species street contains only the streets a car can drive on
   //this variable should be filled with the maximus speed a car is allowed to go on the street
   float maxspeedT;
   float maxspeedB;
   bool in_endgame <- false;
   rgb color <- #blue ;
   float street_width <- 2;
   bool arrow_green <- false;
   bool arrow_red <- false;
 
  
   
   aspect base {   	
   	draw shape color: color width: street_width;   	
   	 // Call draw_arrow_along within the aspect to visualize the direction
    if (maxspeedT = -1) {
        // Only backward direction allowed
        do draw_arrow_along(shape.points, #red);
    } else if (maxspeedB = -1) {
        // Only forward direction allowed
       do draw_arrow_along(shape.points, #green);
    } 
   }
   
action setup_graph {
    // Check if the road is a one-way street
    if (maxspeedT = -1) {
        // If only backward direction is allowed (one-way street)
        color <- #red;
        
        // Reverse the geometry of the road to match the backward direction
        shape <- polyline(reverse(shape.points));
        
    } else if (maxspeedB = -1) {
        // If only forward direction is allowed (one-way street)
        color <- #green;
        
        // No need to reverse geometry, as it aligns with the forward direction
    } else if (maxspeedT > -1 and maxspeedB > -1) {
        // If the street is bidirectional, create a duplicate street in reverse direction
        create street {
            shape <- polyline(reverse(myself.shape.points));
            color <- myself.color;        
          	}
          	
          	}
            
        }
 action draw_arrow_along(list<point> points, rgb color) {
    float arrow_spacing <- 5.0;  // Distance between arrows along the line
    float arrow_size <- 3;     // Size of the arrowhead
    float wing_offset <- 3;     //  Offset to pull back the wings of the arrow
    
    // Traverse the points to place arrows at intervals
    float accumulated_distance <- 0.0;
    loop i from: 0 to: (length(points) - 2) {
        point start <- points[i];
        point end <- points[i + 1];
        
        // Calculate the distance between points
        float segment_length <- distance_to(start, end);
        accumulated_distance <- accumulated_distance + segment_length;
        
        // Draw an arrow if we've reached the arrow spacing threshold
        if accumulated_distance >= arrow_spacing {
            accumulated_distance <- 0.0;  // Reset distance accumulator
            
            // Draw main line segment of the arrow
            //draw line([start, end]) color: color;
            
            // Calculate the direction vector
            point direction <- normalize_vector(start, end);
            
            // Calculate the arrowhead points
            point arrow_tip <- [end.x - 0.5 * direction.x, end.y - 0.5 * direction.y];
           point left_wing <- [arrow_tip.x - (wing_offset + arrow_size) * direction.x + arrow_size * (-direction.y), 
                                arrow_tip.y - (wing_offset + arrow_size) * direction.y + arrow_size * direction.x];
            point right_wing <- [arrow_tip.x - (wing_offset + arrow_size) * direction.x - arrow_size * (-direction.y), 
                                 arrow_tip.y - (wing_offset + arrow_size) * direction.y - arrow_size * direction.x];
            // Draw the arrowhead wings
            draw line([arrow_tip, left_wing]) color: color;
            draw line([arrow_tip, right_wing]) color: color;
        }
    }
}

// Helper function to normalize a vector between two points
action normalize_vector(point p1, point p2) {
    float dx <- p2.x - p1.x;
    float dy <- p2.y - p1.y;
    float length <- sqrt(dx^2 + dy^2);
    return [dx / length, dy / length];
}   

   


// Funktion zur Berechnung und Normalisierung eines Vektors zwischen zwei Punkten
action normalize_vector(point p1, point p2) {
    // Berechne die Differenz (Vektor) zwischen den beiden Punkten
	float dx <- p2.x - p1.x;
	float dy <- p2.y - p1.y;
    // Berechne die Länge des Vektors
	float length <- sqrt(dx^2 + dy^2);
	point normalized_vector <- [0, 0]; // Erstelle einen leeren Punkt als Vektor
	if (length > 0) {
        // Normalisiere den Vektor
		normalized_vector <- [dx / length, dy / length]; // Normierter Vektor
	}
	return normalized_vector;  // Gib den normalisierten Vektor zurück
}
   
   
}

species adress {
   //number of households at adress
   int households;
   //number of parcels received
   int recparcels;
   rgb color <- #red;
   float min_distance_to_adress <- 1e9;
   point closest_street_location;
   point first_vertex;
   point last_vertex;
   graph street_graph;
   point temp_closest_point_on_foot;
   point temp_closest_point_on_driveable_street;
   point closest_point_on_driveable_street;
   point closest_point_on_foot;
   float min_point_distance <- 1e9;
   float min_driveable_point_distance <-1e9;
   bool considered_on_way_to_vehicle <- false; 
   list<street> streets_for_endgame;
   point closest_point_on_walkable_path;
   point closest_point_on_foot_to_point_on_driveable_street;

   aspect base {
        draw circle(2) color: #red;
        draw circle(2) at: closest_point_on_foot color: #green; 
        draw circle(2) at: temp_closest_point_on_driveable_street color: #purple; 
        
        
   } 

//action to find the point on the closest street which ist closest to the adress - resembles a garden gate
   action find_closest_point_on_foot {
   	    street_to_find_entrance nearest_street <- nil;
        float min_street_distance <- 1e9;
        list<point> points_to_choose <-nil; 
        list<street_to_find_entrance> excluded_streets <-nil;       
        //repeated 4 times - the 4 closest streets are considered              
        loop times: 10{              
       ask street_to_find_entrance {
       	   if (not (self in excluded_streets)){
       		
            	float distance_to_adress <- self.location distance_to myself.location;  
            	if (distance_to_adress < min_street_distance) {
                	min_street_distance <- distance_to_adress;
                	nearest_street <- self;
                geometry street_geometry <- self.shape;
                list<point> street_geometry_points_list <- street_geometry.points;
                myself.first_vertex <- street_geometry.points[0].location;
                myself.last_vertex <- street_geometry.points[length(street_geometry_points_list) - 1].location;
                myself.street_graph <- as_edge_graph(nearest_street);
                }
            
            }
        }
        
        if (nearest_street != nil) {
     	  excluded_streets <- excluded_streets + nearest_street;
     	    } 	

     	//create a not visible walker which creates points along a street (1m steps) - the closest point to the adress will be saved 
     	create walker number: 1 {
     		location <- myself.first_vertex.location;
     		}
     		
     	ask walker {
     		bool stillgoing <- true;
     		loop while: stillgoing {
     			do goto target:myself.last_vertex on: myself.street_graph speed: 1.0;
     			float distance_adress <- location distance_to myself.location;
     			
     			if (distance_adress < myself.min_point_distance) {
            	    myself.min_point_distance <- distance_adress;
            	    myself.temp_closest_point_on_foot <- self.location;
           	     }
     			if location = myself.last_vertex.location {
     				stillgoing <- false;
     				points_to_choose <- points_to_choose + myself.temp_closest_point_on_foot;
     				}
     			}
     		do die;
     		}
     	      	       	if self.name = "adress283"{

     	}
     	
     	   min_point_distance <- 1e9;	
     	   min_street_distance <- 1e9;
     	   
 
     	}
     	
     	
     	
	//choosing the closest of all 4 closest points
	closest_point_on_foot <- with_min_of(points_to_choose, each distance_to location);

	}
	
   action find_closest_point_on_driveable_street {
   	    street nearest_street <- nil;
        float min_street_distance <- 1e9;
        list<point> points_to_choose <-nil; 
        list<street> excluded_streets <-nil;       
        //repeated 4 times - the 4 closest streets are considered              
        loop times: 10{              
       ask street {
       	   if (not (self in excluded_streets)){
       		
            	float distance_to_adress <- self.location distance_to myself.location;  
            	if (distance_to_adress < min_street_distance) {
                	min_street_distance <- distance_to_adress;
                	nearest_street <- self;
                geometry street_geometry <- self.shape;
                list<point> street_geometry_points_list <- street_geometry.points;
                myself.first_vertex <- street_geometry.points[0].location;
                myself.last_vertex <- street_geometry.points[length(street_geometry_points_list) - 1].location;
                myself.street_graph <- as_edge_graph(nearest_street);
                }
            
            }
        }
        
        if (nearest_street != nil) {
     	  excluded_streets <- excluded_streets + nearest_street;
     	    } 	
     	//create a not visible walker which creates points along a street (1m steps) - the closest point to the adress will be saved 
     	create walker number: 1 {
     		location <- myself.first_vertex.location;
     		}
     		
     	ask walker {
     		bool stillgoing <- true;
     		loop while: stillgoing {
     			do goto target:myself.last_vertex on: myself.street_graph speed: 1.0;
     			float distance_adress <- location distance_to myself.location;
     			
     			if (distance_adress < myself.min_driveable_point_distance) {
            	    myself.min_driveable_point_distance <- distance_adress;
            	    myself.temp_closest_point_on_driveable_street <- self.location;
           	     }
     			if location = myself.last_vertex.location {
     				stillgoing <- false;
     				points_to_choose <- points_to_choose + myself.temp_closest_point_on_driveable_street;
     				}
     			}
     		do die;
     		}
     	 
     	
     	   min_point_distance <- 1e9;	
     	   min_street_distance <- 1e9;
     	   
     	   
     	}
     	
     	
     	
	//choosing the closest of all 4 closest points
	closest_point_on_driveable_street <- with_min_of(points_to_choose, each distance_to location);

	}	
	
}

species walker skills:[moving] {
            //this species serves only the purpose of finding closest points on lines
    }


species deliverer skills:[moving] {
	float dspeed <- movespeed;
	int maxparcels <- maximum_parcel_carry_capacity;
	int parcels <- maxparcels;
	bool incar <- false;
	rgb color <- #black;
	point the_target <- nil ;
	point the_target_entrance <- nil;
	string the_target_name <- nil;
	float distance_to_next_adress <- nil;
	int needed_parcels;
	bool in_entrance <- false;
	point exitpoint;
	bool on_wayout <- false;
	bool loopstopper <- false;
	path show_me_the_path;
	list<string> already_considered <- [];
	bool after_delivery <- false;
	
	
	aspect base{
		draw file("D:/Master/Master-Thesis/Programme/GAMA/deliverer.svg") at: {location.x, location.y, 10} size: 10 color: color rotate:180;
		//draw circle(100) at:the_target_entrance color: #orange; 
		draw show_me_the_path color: #yellow width: 5; 
	} 
	//when in a vehicle the deliverer drives 
	reflex drive_on when: incar {
		ask vehicle {
			if length(adress) !=0{
			do drive;			
			}
		}
	}
	
	 // initial reflex to find new target
    reflex initial_target when: the_target = nil and not incar and not on_wayout{
        do find_nearest_adress;
    }
	// Method to set the nearest address(actually the nearest enter point!!) as destination
	action find_nearest_adress {
		
		if length(adress) != 0 and the_target = nil{
	
        adress nearest_adress <- tsp_nodes[optimal_path[0]];        
        float min_distance <- 1e9;  // setting a very big number as min_distance
        
//        ask adress {
//           //alt float distance_to_adress <- closest_point_on_foot.location distance_to (myself.location);
//           path path_to_adress <- path_between(the_walking_graph, myself.location, closest_point_on_foot.location);
//			
//			if path_to_adress != nil {
//				
//			    float total_length <- 0.0;
// 			   list<geometry> segments <- path_to_adress.segments;
//
// 			   loop segment over: segments {
// 			       float segment_length <- segment.perimeter;
// 			       total_length <- total_length + segment_length;
//  			  }
//
//  			  float distance_to_adress <- total_length;
//
// 			   if (distance_to_adress <= min_distance) {
// 			       min_distance <- distance_to_adress;
// 			       nearest_adress <- self;
//  			      myself.needed_parcels <- self.households;
// 			   }
//			}
//			
//			
//			}
			
		
        if (nearest_adress != nil) {
            the_target <- nearest_adress.location; 
            the_target_name <- nearest_adress.name;
            the_target_entrance <- nearest_adress.closest_point_on_foot.location; 
			
        
        
        
        
        show_me_the_path <- path_between(the_walking_graph, location, nearest_adress.closest_point_on_foot.location);
        float distance_to_next_adress_on_path <- 0;
		list<geometry> segments <- show_me_the_path.segments;
		loop segment over: segments {
            float segment_length <- segment.perimeter;  // Berechne die Länge des Segments
            distance_to_next_adress_on_path <- distance_to_next_adress_on_path + segment_length;  // Füge die Segmentlänge zur Gesamtlänge hinzu
        }
        distance_to_next_adress <- distance_to_next_adress_on_path;
        if endgame {
        	distance_until_driving <-distance_to_next_adress +1;
       
        }
}else{
	write "Error: no nearest adress!";
	do return_to_vehicle;
}
        
        
        }
	}
	//this action finds the nearest adress and prevents the changing of the target -->loopstopper
	action find_nearest_adress_during_return {
		if length(adress) != 0 and the_target = nil{
        adress nearest_adress <-tsp_nodes[optimal_path[0]];
        float min_distance <- 1e9;  // setting a very big number as min_distance
        
//        ask adress {
//           //alt float distance_to_adress <- closest_point_on_foot.location distance_to (myself.location);
//           path path_to_adress <- path_between(the_walking_graph, myself.location, closest_point_on_foot);
//			
//			if path_to_adress != nil {
//			    float total_length <- 0.0;
// 			   list<geometry> segments <- path_to_adress.segments;
//
// 			   loop segment over: segments {
// 			       float segment_length <- segment.perimeter;
// 			       total_length <- total_length + segment_length;
//  			  }
//
//  			  float distance_to_adress <- total_length;
//
// 			   if (distance_to_adress <= min_distance) {
// 			       min_distance <- distance_to_adress;
// 			       nearest_adress <- self;
//  			      myself.needed_parcels <- self.households;
// 			   }
//			}
//			
//			
//			}
        if (nearest_adress != nil) {
            the_target <- nearest_adress.location; 
            the_target_name <- nearest_adress.name;
            the_target_entrance <- nearest_adress.closest_point_on_foot.location; 
        
        
 		path path_to_adress <- path_between(the_walking_graph, location, nearest_adress.closest_point_on_foot.location);
           float distance_to_next_adress_on_path <- 0;
		list<geometry> segments <- path_to_adress.segments;
		loop segment over: segments {
            float segment_length <- segment.perimeter;  // Berechne die Länge des Segments
            distance_to_next_adress_on_path <- distance_to_next_adress_on_path + segment_length;  // Füge die Segmentlänge zur Gesamtlänge hinzu
        }
        distance_to_next_adress <- distance_to_next_adress_on_path;
        if endgame {
        	if not after_delivery{
        	distance_until_driving <-distance_to_next_adress +1;
        	
        	}
        }
        if needed_parcels <= parcels {
        	loopstopper <-true;
        }
        }else{
        	write "Error: no nearest adress during return!";
        	}
	}
	}
	
	
	//considering return to vehicle
	//return to vehicle if the deliverer has not enough parcels or the next adress is too far away - the deliverer can not return to the vehicle until he is back at the enter point or walk
	reflex considering when: incar = false {
		        	
		if ((needed_parcels > parcels or distance_to_next_adress > distance_until_driving) and on_wayout=false){
			the_target <- nil; 
            the_target_name <- nil;
            the_target_entrance <- nil; 
			do return_to_vehicle;
		}
		else {
			do walk;
		}
		
	}
	
	action return_to_vehicle {
		// calculating distance to vehicle
        float distance_to_vehicle <- location distance_to any (vehicle).vehicle_substitute.location;
        point vehicle_location <- location of one_of (vehicle).vehicle_substitute;
        point true_vehicle_location <- location of one_of (vehicle);
       
 if (distance_to_vehicle > 0) {
 	
        do goto target: any(vehicle).vehicle_substitute.location on: the_walking_graph speed: dspeed;
        
    // weil der deliverer ansonsten bei der letzten adresse am walking graph stehen bleibt und goto einfach nicht ausführt.
    if length(adress) = 0{    
    do goto target: any(vehicle).vehicle_substitute.location on: the_walking_graph speed: dspeed; 
    total_time <- total_time +1;
    }   
    //}

    
}

        	//find a new target as there might be an opportunity on the way to the vehicle
        	if length (adress) !=0 {
        	do find_nearest_adress_during_return;
        	
        	}

        	// checking if deliverer is at vehicle and refill the parcels -maybe there should be time passing-
        	if (location = vehicle_location) {
            	parcels <- maxparcels;
            	
            	do find_nearest_adress;
            	if length(adress) = 0{
            		location <- true_vehicle_location;
            		incar<-true;
            	}
        	}
        	//if the deliverer enters the car time will pass, deliverer enters if next adress is too far away to walk
        	if (distance_to_next_adress >= distance_until_driving and location = vehicle_location){
        		incar <-true;
        		location <- true_vehicle_location;
        		loopstopper <- false;  
        		total_time <- total_time + enter_exit_time;
        		time_for_vehicle_interaction <- time_for_vehicle_interaction+ enter_exit_time;       		 						
        	}
        	
	}
	 // movement to the target on the walkable streets
    action walk {
    	if (distance_to_next_adress < distance_until_driving and parcels >= needed_parcels and  in_entrance = false) {
    		
    		do goto target:the_target_entrance on: the_walking_graph speed:dspeed;
    		
        
        // checking if deliverer is at target location and entering the entrance
        if location = the_target_entrance {
        	        	in_entrance <- true;
        	        	
        }
        }
        //moving without a graph in a straigth line
        if in_entrance = true and on_wayout = false{
        	do goto target: the_target speed: dspeed;
        	
        }
        //returning to the entrance after delivery
        if in_entrance = true and on_wayout = true {      	
        	do goto target: exitpoint speed: dspeed;
        	//resetting 
        	if location = exitpoint {
        		exitpoint <- nil;
        		on_wayout <- false;
        		in_entrance <- false;
        		after_delivery <- true;
//        		if length(adress) != 0{
//	        	do find_nearest_adress;
//   		     	if distance_to_next_adress > distance_until_driving {
// 		       	distance_until_driving <-distance_until_driving_fixed;
// 		       	} 		       	
// 		       	}
 		       	after_delivery <- false;
 		       	the_target <- nil;
   		     	the_target_name <- nil;
        		do find_nearest_adress;
        	}
        }
        if (location = the_target) {
        	// deliver parcels if at adress
        	do deliver_parcel;  
        }
    }
    
    action walk_without_checking {
    	//should be used when going to a adress without searching for better (closer adresses) to prevent endles loops
    	if (distance_to_next_adress < distance_until_driving and parcels >= needed_parcels and  in_entrance = false) {
    		
    		do goto target:the_target_entrance on: the_walking_graph speed:dspeed;
            if location = the_target_entrance {
        	        	in_entrance <- true;
        	        	
        }
        }
        if in_entrance = true and on_wayout = false{
        	do goto target: the_target speed: dspeed;
        	
        }
        if in_entrance = true and on_wayout = true {
        	do goto target: exitpoint speed: dspeed;
        	if location = exitpoint {
        		exitpoint <- nil;
        		on_wayout <- false;
        		in_entrance <- false;
        		after_delivery <- true;
//        		if length(adress) != 0{
//	        	do find_nearest_adress;
//   		     	if distance_to_next_adress > distance_until_driving {
// 		       	distance_until_driving <-distance_until_driving_fixed;
// 		       	} 		       	
// 		       	}
 		       	after_delivery <- false;
 		       	the_target <- nil;
   		     	the_target_name <- nil;
    	        if endgame {        	        	
	        	endgame <-false;
       			}
        		}
        		
        }
        if (location = the_target) {
        	
            do deliver_parcel;  // deliver parcels if destination is distance_to_next_adressed
        }
    }
    
    // deliver parcels and delete adress 
    action deliver_parcel {
    	// set exitpoint so the deliverer can return to the street
    	exitpoint <- the_target_entrance;
    	on_wayout <-true;
        ask adress {
        	if (self.location = myself.location){
            myself.parcels <- myself.parcels - households;
            recparcels <- recparcels + households;
            //looking inside the time matrix on position of the number of households (-1 because 1 household ist at line 0 )
            string stringfrommatrix <- time_matrix[1,households-1];
            float floatfrommatrix <- float(stringfrommatrix);
            //time is in min --> conversion to s
            float time_per_address <- floatfrommatrix * 60;
            
            time_for_delivery <- time_for_delivery + time_per_address;
            //adding the time
            total_time <- total_time + time_per_address;
			// deleting of adress after delivery
            if (recparcels >= households) {
                do die;
            }
            
            }
        }
        optimal_path <- optimal_path - optimal_path[0];
        write optimal_path;
        // find next adress
        the_target <- nil;
        the_target_name <- nil;
        if endgame {        	        	
        	endgame <-false;
        }
       }
    //gotocar serves the purpose of positioning the deliverer at the car on the street. the time walking there will be added.
   action gotocar {
   	bool ongoing <- true;
   	loop while: ongoing {
   		float move_distance <- min(location distance_to any(vehicle), dspeed * step);
   	    point carlocation <- any (vehicle);
     point move_step <- (location + ((carlocation - location)*(move_distance / (location distance_to carlocation))));
       	        
        //moving the agent to the next location
        location <- move_step;
        //total_time <- total_time +1;
        if location = carlocation {
        	ongoing <- false;
        	incar <-true;
        	//total_time <- total_time + enter_exit_time;
        	//time_for_vehicle_interaction <- time_for_vehicle_interaction+ enter_exit_time;
        }
         
   }
   
   }
   
  
   
}






species vehicle skills: [moving] {
	float vspeed <- vehiclespeed;
	rgb color <- #orange;
	point next_parking_spot <- nil;
	float vehicle_Start_min_dist <- 1e9;
	float vehicle_parking_spot_min_dist <- 1e9;
	point closest_street_location;
    point first_vertex_vehicle;
    point last_vertex_vehicle;
    graph street_graph_vehicle; 
    float min_start_point_distance <- 1e9;
	point temp_closest_start_point_on_street;
	point closest_starting_point_on_street;
	list<point> possible_parking_spots;
	list<geometry> exact_possible_parking_spots;
	list<point> choosen_fallback_parking_spots;
	bool is_forward <- true;
	point previous_location <- nil;
	bool first_parking_spot_reached <- false;
	list<street>streets_for_endgame;
	path path1;
	path path2;
	string next_parking_spot_name;
	// Temporary point on the walking graph to represent the vehicle's location
    point vehicle_substitute <- nil;  
    
	
	
	aspect base{
		draw file("D:/Master/Master-Thesis/Programme/GAMA/car.svg") size: 10 color: color rotate:180; 
		if (vehicle_substitute != nil) {
            draw circle(5) at: vehicle_substitute color: #blue;  // Visualize the substitute point
        }
		if (next_parking_spot != nil) {
            draw circle(20) at: {next_parking_spot.x,next_parking_spot.y, 9} color: #blue border: #black;  // Darstellung des nächsten Parkplatzes
            loop i over:choosen_fallback_parking_spots{
            draw circle (5) at: i color: #green;
            
            }
        draw path1.shape color: #green width: 5;
        draw path2.shape color: #purple width: 5;
        }
	} 
	
	
	//finding the nearest road to the starting point and positioning the vehicle there
	action find_start {
		street nearest_starting_street <- nil;
        float min_starting_street_distance <- 1e9;
        list<point> points_to_choose_start; 
        list<street> excluded_starting_streets;       
        //looping 4 times to consider the 4 nearest roads               
        loop times: 4{              
       ask street {
       	   if (not (self in excluded_starting_streets)){
       		
            	float distance_to_start <- self.location distance_to any(start_point);  
            	if (distance_to_start < min_starting_street_distance) {
                	min_starting_street_distance <- distance_to_start;
                	nearest_starting_street <- self;
                geometry street_geometry_vehicle <- self.shape;
                list<point> street_geometry_vehicle_points_list <- street_geometry_vehicle.points;
                myself.first_vertex_vehicle <- street_geometry_vehicle.points[0].location;
                myself.last_vertex_vehicle <- street_geometry_vehicle.points[length(street_geometry_vehicle_points_list) - 1].location;
                myself.street_graph_vehicle <- as_edge_graph(nearest_starting_street);
                }
            
            }
        }
        
        if (nearest_starting_street != nil) {
     	  excluded_starting_streets <- excluded_starting_streets + nearest_starting_street;
     	    } 	
     	 
     	create walker number: 1 {
     		location <- myself.first_vertex_vehicle.location;
     		}
     		
     	ask walker {
     		bool stillgoing <- true;
     		loop while: stillgoing {
     			do goto target:myself.last_vertex_vehicle on: myself.street_graph_vehicle speed: 1.0;
     			float distance_start <- location distance_to any(start_point);
     			
     			if (distance_start < myself.min_start_point_distance) {
            	    myself.min_start_point_distance <- distance_start;
            	    myself.temp_closest_start_point_on_street <- self.location;
           	     }
     			if location = myself.last_vertex_vehicle.location {
     				stillgoing <- false;
     				points_to_choose_start <- points_to_choose_start + myself.temp_closest_start_point_on_street;
     				}
     			}
     		do die;
     		}
     	 
     	
     	   min_start_point_distance <- 1e9;	
     	   min_starting_street_distance <- 1e9;
     	        	   
     	}
     	    	
     	
closest_starting_point_on_street <- with_min_of(points_to_choose_start, each distance_to any(start_point));
location <- closest_starting_point_on_street.location;
}



// Action to find parking spots
action find_parking_spots {
    possible_parking_spots <- nil;
    exact_possible_parking_spots <- nil;
    vehicle_parking_spot_min_dist <- 1e9;
    write optimal_path[0];
    write tsp_nodes[optimal_path[0]];
    write tsp_nodes[optimal_path[0]].closest_point_on_driveable_street;
	next_parking_spot <- tsp_nodes[optimal_path[0]].closest_point_on_driveable_street;
//     Start at a random parking spot if first parking spot has not been reached
//    if (first_parking_spot_reached = false) {    
//        next_parking_spot <- adress[simulation_counter].closest_point_on_driveable_street;   
//        
//        path path_to_next_spot  <- path_between(the_graph, location, next_parking_spot);
//     	
//     	if (path_to_next_spot = nil) {
//            write "No valid path to the first parking spot. Nearest pointwill be chosen.";
//            // Find the nearest connected non-oneway street
//            street nearest_connected_street <- nil;
//            float min_distance <- 1e9;
//			list<street> nearest_streets_ini <- [];
//			
//			ask adress[simulation_counter] {               
//                float min_distance;
//                loop times: 4 {
//                    float nearest_distance <- 1e9;
//                    street nearest_street <- nil;
//                    ask street {
//                        float distance_to_adress <- location distance_to myself.closest_point_on_driveable_street;
//                        if (distance_to_adress < nearest_distance and not(self in nearest_streets_ini)) {
//                            nearest_distance <- distance_to_adress;
//                            nearest_street <- self;
//                        }
//                    }
//                    if nearest_street != nil {
//                        nearest_streets_ini <- nearest_streets_ini + nearest_street;
//                    }
//                
//               }
//                           
//                }   
//			
//			float end_distance <- 1e9; // Initialize with a large value
//            // Iterate over all streets to find a suitable connected endpoint
//            loop i over: nearest_streets_ini {
//                if (i.maxspeedT > -1 and i.maxspeedB > -1) { // Exclude one-way streets in the wrong direction
//                    list<point> street_points <- i.shape.points;
//
//                    // Check the start and end points of the street
//                    point start_point_normal_street <- street_points[0];
//                    point end_point_normal_street <- street_points[length(street_points) - 1];
//
//                    // Determine the connection distance
// 					// Calculate the path distance to the start point
//					path path_to_start <- path_between(the_walking_graph, next_parking_spot, start_point_normal_street);
//                    float start_distance <- 1e9; // Initialize with a large value
//                    if (path_to_start != nil) {
//                        start_distance <- 0.0;
//                        list<geometry> segments <- path_to_start.segments;
//                        loop segment over: segments {
//                            start_distance <- start_distance + segment.perimeter;
//                        }
//                    }
//
//                    // Calculate the path distance to the end point
//                    path path_to_end <- path_between(the_walking_graph, next_parking_spot, end_point_normal_street);
//                    
//                    if (path_to_end != nil) {
//                        end_distance <- 0.0;
//                        list<geometry> segments <- path_to_end.segments;
//                        loop segment over: segments {
//                            end_distance <- end_distance + segment.perimeter;
//                        }
//                    }
//                    
//                    if (start_distance < min_distance) {
//                        min_distance <- start_distance;
//                        nearest_connected_street <- self;
//                        next_parking_spot <- start_point_normal_street;
//                    }
//
//                    if (end_distance < min_distance) {
//                        min_distance <- end_distance;
//                        nearest_connected_street <- self;
//                        next_parking_spot <- end_point_normal_street;
//                    }
//                }
//            }
//        }
//        
//     	         
//		if (next_parking_spot in point_to_address_name) {
//        next_parking_spot_name <- point_to_address_name[next_parking_spot];
//  		  } else {
//        next_parking_spot_name <- "Unknown"; 
//        
//    }   
//    
//      
//    } else {
//        // Find possible parking spots near addresses within walking distance
//       ask adress {
//            path distance_to_enter <- path_between(the_walking_graph, self.closest_point_on_driveable_street, self.closest_point_on_foot);
//            float total_length <- 0.0;
//            if distance_to_enter != nil{
//            list<geometry> segments <- distance_to_enter.segments;
//            loop segment over: segments {
//                total_length <- total_length + segment.perimeter;
//            }
//
//            if total_length < distance_until_driving {
//                myself.possible_parking_spots <- myself.possible_parking_spots + closest_point_on_driveable_street;                        
//            }
//            }
//         
//        }
//
//        // Calculate the closest parking spot based on path length
//        if possible_parking_spots != [] {
//            loop i over: possible_parking_spots {
//                path path_to_next_spot <- path_between(the_graph, location, i);
//if (path_to_next_spot = nil) {
//	possible_parking_spots <- possible_parking_spots - i;
//}
//                if (path_to_next_spot != nil) {
//                    float total_length <- 0.0;
//                    list<geometry> segments <- path_to_next_spot.segments;
//
//                    // Add the length of all segments except the last one
//                    loop segment over: segments {
//                        if segment != segments[length(segments) - 1] {
//                            total_length <- total_length + segment.perimeter;
//                        }
//                    }
//
//                    // Handle the last segment only up to the target point
//			        if length(segments) > 0 {
//  			          geometry last_segment <- segments[length(segments) - 1];
//  			          list<point> last_segment_points <- last_segment.points;
//   			         point start_node_of_last_segment <- last_segment_points[0];
//  			          float segment_distance_to_target <- start_node_of_last_segment distance_to i;
//
//    			        // Add the distance from the start node of the last segment to the target point
//    			        total_length <- total_length + segment_distance_to_target;
//    			    }
//
//                    // Update the next parking spot if this one is closer
//                    if (total_length < vehicle_parking_spot_min_dist) {
//                        vehicle_parking_spot_min_dist <- total_length;
//                        next_parking_spot <- i;
//                        next_parking_spot_name <- point_to_address_name[i];
//                    }
//                }
//            }
//           
//        }
// 
//		if possible_parking_spots = [] {
//		
//			endgame <-true;
//            // If no possible parking spots found, look for fallback spots
//            list<point> fallback_parking_spots <- nil;
//            list<point> fallback_parking_spots_before <- nil;
//            // Find the 4 nearest streets to each closest_point_on_foot
//            ask adress {
//                list<street> nearest_streets <- [];
//                float min_distance;
//                loop times: 4 {
//                    float nearest_distance <- 1e9;
//                    street nearest_street <- nil;
//                    ask street {
//                        float distance_to_adress <- location distance_to myself.closest_point_on_foot;
//                        if (distance_to_adress < nearest_distance and not(self in nearest_streets)) {
//                            nearest_distance <- distance_to_adress;
//                            nearest_street <- self;
//                        }
//                    }
//                    if nearest_street != nil {
//                        nearest_streets <- nearest_streets + nearest_street;
//                    }
//                
//                           
//                }     
//                // Find the closest points on these streets to the closest_point_on_foot
//                loop each_street over: nearest_streets {
//                    float min_fallback_distance <- 1e9;
//                    point best_fallback_spot <- nil;
//                    list<point> street_points <- each_street.shape.points;
//                    
//                   loop street_point over: street_points {
//                        path path_to_foot <- path_between(the_walking_graph, street_point, closest_point_on_foot);
//                        if path_to_foot != nil {
//                            float path_length <- 0.0;
//                            list<geometry> segments <- path_to_foot.segments;
//                            loop segment over: segments {
//                                path_length <- path_length + segment.perimeter;
//                            }
//                             if path_length < min_fallback_distance {
//                                min_fallback_distance <- path_length;
//                                best_fallback_spot <- street_point;
//                            }                                                       
//                            }
//                            }
//                            if best_fallback_spot != nil {
//                        fallback_parking_spots_before <- fallback_parking_spots_before + best_fallback_spot;
//                    }
//                            }
//                            
//                   loop each_pos_point over: fallback_parking_spots_before {
//                   	point best_fallback_spot <- nil;
//                    float min_fallback_distance <- 1e9;
//                     path path_to_foot <- path_between(the_walking_graph, each_pos_point, closest_point_on_foot);
//                        if path_to_foot != nil {
//                            float path_length <- 0.0;
//                            list<geometry> segments <- path_to_foot.segments;
//                            loop segment over: segments {
//                                path_length <- path_length + segment.perimeter;
//                            }
//                             if path_length < min_fallback_distance {
//                                min_fallback_distance <- path_length;
//                                best_fallback_spot <- each_pos_point;
//                            }                                                       
//                            }                         
//                            if best_fallback_spot != nil {
//                        fallback_parking_spots <- fallback_parking_spots + best_fallback_spot;
//                    }
//                    }
//                    }
//                            
//                    // Select the fallback parking spot that is closest to the vehicle's current location
//            float min_distance_to_vehicle <- 1e9;
//            point chosen_fallback_spot <- nil;
//                    
//                    loop fallback_spot over: fallback_parking_spots {
//                    path path_to_vehicle <- path_between(the_graph, location, fallback_spot);
//                    
//                     if path_to_vehicle != nil {
//                    float total_path_length <- 0.0;
//                    list<geometry> segments <- path_to_vehicle.segments;
//					
//					if segments != []{	
//                    // Calculate the length of the path
//                    loop segment over: segments {
//                    	if segment != segments[length(segments) - 1] {
//                        total_path_length <- total_path_length + segment.perimeter;
//                        }
//                    }
//                    
//                     // Now handle the last segment only up to the target point
//                    geometry last_segment <- segments[length(segments) - 1];
//                    list<point> last_segment_points <- last_segment.points;
//                    point start_node_of_last_segment <- last_segment_points[0];
//                    float segment_distance_to_target <- start_node_of_last_segment distance_to fallback_spot;
//                    
//                     // Add the distance from the start node of the last segment to the target point
//                    total_path_length <- total_path_length + segment_distance_to_target;
//
//                    // Check if this path is the shortest path to a fallback spot
//                    if total_path_length < min_distance_to_vehicle {
//                        min_distance_to_vehicle <- total_path_length;
//                        chosen_fallback_spot <- fallback_spot;
//                    }
//                    }
//                }
//                    }
//                    
//              if chosen_fallback_spot != nil {
//              	
//                next_parking_spot <- chosen_fallback_spot;
//                if (next_parking_spot in point_to_address_name) {
//        			next_parking_spot_name <- point_to_address_name[next_parking_spot];
//  		  			} else {
//       				next_parking_spot_name <- "Unknown3";       				
//       				}
//            } else {
//                // If no fallback spot found, choose a random address's closest point as a last resort
//                
//                    next_parking_spot <- any(adress).closest_point_on_driveable_street;
//                    if (next_parking_spot in point_to_address_name) {
//        			next_parking_spot_name <- point_to_address_name[next_parking_spot];
//  		  			} else {
//       				next_parking_spot_name <- "Unknown4";       				
//       				}
//                
//            }
//                    }
//                                  
//                 }
                 
                 }                                         	
		
			
	
	action drive  {
		float nearest_maxspeed <- 0 #km / #h;
		float min_distance <- 1e9;
        street nearest_street <- nil;
        if next_parking_spot = nil{
        	do find_parking_spots;
          	}
        	
        // find nearest street
        ask street {
            float distance_to_vehicle <- self.location distance_to (location);
            if (distance_to_vehicle < min_distance) {
                min_distance <- distance_to_vehicle;
                nearest_street <- self;
            }
        }

        if (nearest_street != nil) {
            nearest_maxspeed <- nearest_street.maxspeedT;
        }
		 vspeed <- nearest_maxspeed #km / #h;
		 
		do goto(target: next_parking_spot, on:the_graph, speed: vspeed);
	ask deliverer {
		location <- myself.location;
	}
		if (next_parking_spot = self.location ) {
			if first_parking_spot_reached =false{
				time_to_area <- total_time;
				
			}
			first_parking_spot_reached <- true;
			next_parking_spot <- nil ;
			// Create a substitute point on the walking graph near the vehicle's location
        	list<point> Possible_points;
        	list<walkablestreet> excluded_walkable_streets ;
        	float min_walkable_street_distance <- 1e9;
        	walkablestreet nearest_walkable_street;
       loop times:4{ 	
        ask walkablestreet {
       	   if (not (self in excluded_walkable_streets)){
       		
            	float distance_to_vehicle <- self.location distance_to myself.location;  
            	if (distance_to_vehicle < min_walkable_street_distance) {
                	min_walkable_street_distance <- distance_to_vehicle;
                	nearest_walkable_street <- self;
                }
            
            }              
        }
        
        if (nearest_walkable_street != nil) {
     	  excluded_walkable_streets <- excluded_walkable_streets + nearest_walkable_street;
     	    }
     	    
     	    }
     	loop i over: excluded_walkable_streets {
     		 	
     	        geometry street_geometry <- i.shape;
                list<point> street_geometry_points_list <- street_geometry.points;
                point first_vertex <- street_geometry.points[0].location;
                point last_vertex <- street_geometry.points[length(street_geometry_points_list) - 1].location;
                graph street_graph <- as_edge_graph(nearest_walkable_street);
     	 
     	//create a not visible walker which creates points along a street (1m steps) - the closest point to the adress will be saved 
     	create walker number: 1 {
            if first_vertex != nil {
                location <- first_vertex; // Set walker to the first vertex
            } else {
                write "Error: First vertex is nil!";
                do die;
            }
        }
     		
     	ask walker {
     		bool stillgoing <- true;
     		float min_point_distance <- 1e9;
     		point temp_closest_point_to_vehicle;
     		loop while: stillgoing {
     			do goto target:last_vertex on:street_graph speed: 1.0;
     			float distance_vehicle <- location distance_to myself.location;
     			
     			if (distance_vehicle < min_point_distance) {
            	    min_point_distance <- distance_vehicle;
            	    temp_closest_point_to_vehicle <- self.location;
           	     }
     			if location = last_vertex.location {
     				stillgoing <- false;
     				Possible_points <- Possible_points + temp_closest_point_to_vehicle;
     				}
     			}
     		do die;
     		}     	      	       	     	
     	}   
     	
     	
     	
     	
     	
			//choosing the closest of all 4 closest points
			vehicle_substitute <- with_min_of(Possible_points, each distance_to location);             	
			ask deliverer{
				location <- myself.vehicle_substitute;
				incar <- false;
				total_time <- total_time + enter_exit_time;
				time_for_vehicle_interaction <- time_for_vehicle_interaction+ enter_exit_time;
				parcels <- maxparcels;
				the_target <- nil;
				the_target_name <- nil;
				
			}
		}
	}
	
	action return_to_base {
		do goto(target: closest_starting_point_on_street, on:the_graph, speed: vspeed);
		time_from_area <- time_from_area + 1;
	ask deliverer {
		location <- myself.location;
		
	}
	}
	
	    action determine_direction {
        // Get the current street the vehicle is on
        float nearest_maxspeed <- 0 #km / #h;
		float min_distance <- 1e9;
        street nearest_street <- nil;
        ask street {
            float distance_to_vehicle <- self.location distance_to (location);
            if (distance_to_vehicle < min_distance) {
                min_distance <- distance_to_vehicle;
                nearest_street <- self;
            }
        }

        street current_street <- nearest_street;
        if (current_street != nil and previous_location = nil) {
        	vspeed <- current_street.maxspeedT;
        }

        if (current_street != nil and previous_location != nil) {
            // Extract the first and last points of the street geometry
            geometry street_geometry <- current_street.shape;
            list<point> street_points <- street_geometry.points;
            point start_vertex <- street_points[0];  // First point of the street
            point end_point <- street_points[length(street_points) - 1];  // Last point of the street

            // Calculate the movement vector of the vehicle
            point movement_vector <- normalize_vector(self.location, previous_location);

            // Calculate the direction vector of the street
            point street_vector <- normalize_vector(end_point, start_vertex);

            // Check if the vehicle is moving in the same direction as the street (forward)
            float dot_product_value <-movement_vector.x *street_vector.x + movement_vector.y *street_vector.y ;
             // Check if the vehicle is moving in the same direction as the street (forward)
            if (dot_product_value > 0) {
                is_forward <- true;
                vspeed <- current_street.maxspeedT;  // Set speed to forward direction speed
            } else {
                is_forward <- false;
                vspeed <- current_street.maxspeedB;  // Set speed to backward direction speed
            }
        }

        // Update the previous location for the next tick
        previous_location <- self.location;
    }
    
    
    action normalize_vector(point p1, point p2) {
    // Berechne die Differenz (Vektor) zwischen den beiden Punkten
	float dx <- p2.x - p1.x;
	float dy <- p2.y - p1.y;
    // Berechne die Länge des Vektors
	float length <- sqrt(dx^2 + dy^2);
	point normalized_vector <- [0, 0]; // Erstelle einen leeren Punkt als Vektor
	if (length > 0) {
        // Normalisiere den Vektor
		normalized_vector <- [dx / length, dy / length]; // Normierter Vektor
	}
	return normalized_vector;  // Gib den normalisierten Vektor zurück
}
	
	
}
species startpunkt {
	float vspeed;
	rgb color <- #green;
	
	aspect base{
		draw shape color: color border: color;
	} 
}


//creating an experiment
experiment Flyer_delivery type: gui {
	//changeable parameters
	parameter "Shapefile for the adresses:" var: adress_shapefile category: "GIS" ;
	parameter "Shapefile for the roads:" var: street_shapefile category: "GIS" ;
	parameter "Shapefile for the stating point:" var: start_point category: "GIS" ;
	parameter "Speed of deliverer:" var: movespeed ;
	parameter "Max parcel capacity" var: maximum_parcel_carry_capacity ;
	parameter "Distance when deliverer desides to drive" var:distance_until_driving_fixed;
	
	output {
		display city_display type: opengl {
			//image file: background_image;
			
			
			species deliverer aspect: base ;
			species building aspect: base;	
			species adress aspect: base ;
			species walkablestreet aspect: base;
			species street aspect: base ;
			species startpunkt aspect: base;	
			species vehicle aspect: base;
			
			
				
		}
		//monitoring some stats for debugging and evaluation
		monitor Deliverer_Status value: (one_of(deliverer).incar ? "In Vehicle" : "On Foot") refresh:every(1 #s);
        monitor Next_Parking_Spot value: (one_of(vehicle).next_parking_spot) refresh:every(1#s);
        monitor position_deliverer value: (one_of(deliverer).location) refresh:every(1#s);
        monitor position_vehicle value: (one_of(vehicle).location) refresh:every(1#s);
        monitor needet_parcels value: (one_of(deliverer).needed_parcels) refresh:every(1#s);
        monitor parcels value: (one_of(deliverer).parcels) refresh:every(1#s);
        monitor next_adress value: (one_of(deliverer).the_target) refresh:every(1#s);
        monitor "Total Time (Seconds)" value: total_time refresh: every(1#s);
        monitor target_adress value: (one_of(deliverer).the_target) refresh:every(1#s);
        monitor target_adress_name value: (one_of(deliverer).the_target_name) refresh:every(1#s); 
        monitor on_wayout value: (one_of(deliverer).on_wayout) refresh:every(1#s);
        monitor Next_Parking_Spot_Adress_Name value: (one_of(vehicle).next_parking_spot_name) refresh:every(1#s);
              
      
      //display chart_display {
      	//chart "Needed Time" type: series {
        //data "Total Time" value: total_time;

        
      //  }
       
   // }
    }
    
}

experiment All_first_parking_spots type: batch repeat: 1 until: fertig keep_seed: true{

   parameter "Erster Parkplatz" var:simulation_counter min: 0 max: 208 step: 1;   
   
   // Declare the results list at the experiment level
    list<float> results <- [];

   output {
   		display city_display type: opengl {
			species deliverer aspect: base ;
			//species building aspect: base;	
			species adress aspect: base ;
			species walkablestreet aspect: base;
			//species street aspect: base ;
			//species startpunkt aspect: base;	
			species vehicle aspect: base;
			}
			
			}

   permanent {
   	
    display Time background: #white refresh: true {
    	chart "Time needed" type: xy {
	        data "Starting point" value: {simulation_counter, total_time};	       
	    }
    }
    display path_length background: #white refresh: true {
    	chart "Length traveled" type: xy {
	        data "Starting point" value: {simulation_counter, total_path_traveled_length};	       
	    }
    }
}
}
experiment Random_first_parking_spots type: batch repeat:1 until: fertig keep_seed: false{
	list<int> random_numbers <- [];
 init {
 	loop times:50{
 		random_numbers <- random_numbers + rnd(0, 722);
 		
 	}
 	write "Generated random numbers: " + random_numbers;
 }   

   parameter "Erster Parkplatz" var:simulation_counter among:random_numbers ;
  // parameter "Erster Parkplatz" var:simulation_counter min:0 max: 209 step:1 ;
   //parameter "Maximale Distanz zu Fuß zur Adresse" var:distance_until_driving_fixed min: 50 max: 500 step: 50;
   

   permanent {
   	
    display Time background: #white refresh: true {
    	chart "Time needed" type: xy {
	        data "Starting point" value: {simulation_counter, total_time};	       
	    }
    }
    display path_length background: #white refresh: true {
    	chart "Length traveled" type: xy {
	        data "Starting point" value: {simulation_counter, total_path_traveled_length};	       
	    }
    }
}

}

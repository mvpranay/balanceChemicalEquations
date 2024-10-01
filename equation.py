import re
import numpy as np
from sympy import Matrix
from pulp import LpProblem, LpVariable, LpInteger, lpSum, LpStatus, value

class Equation:
    def __init__(self, equation_str):
        self.equation_str = equation_str
        # Parse the equation into reactants and products
        self.reactants, self.products = self._parse_equation()
        # Store element counts for reactants and products
        self.reactant_element_counts = [self._count_atoms(compound) for compound in self.reactants]
        self.product_element_counts = [self._count_atoms(compound) for compound in self.products]
        # Create the coefficient matrix
        self.coefficient_matrix = self._create_coefficient_matrix()

    def _parse_equation(self):
        # Split the equation by " = " to separate reactants and products
        reactants_str, products_str = self.equation_str.split(" = ")
        reactants = reactants_str.split(" + ")
        products = products_str.split(" + ")
        return reactants, products

    def _count_atoms(self, compound):
        # Use regex to find elements and their counts in the compound
        matches = re.findall(r'([A-Z][a-z]*)(\d*)', compound)
        element_counts = {}
        for element, count in matches:
            count = int(count) if count else 1
            element_counts[element] = element_counts.get(element, 0) + count
        return element_counts

    def _get_all_elements(self):
        # Get the set of all unique elements in the reactants and products
        elements = set()
        for compound in self.reactants + self.products:
            matches = re.findall(r'([A-Z][a-z]*)', compound)
            elements.update(matches)
        return sorted(list(elements))

    def _create_coefficient_matrix(self):
        # Create a matrix with rows representing elements and columns representing reactants/products
        elements = self._get_all_elements()
        num_compounds = len(self.reactants) + len(self.products)
        matrix = []

        for element in elements:
            row = []
            # Add atom counts for each reactant (positive values)
            for counts in self.reactant_element_counts:
                row.append(counts.get(element, 0))
            # Add atom counts for each product (negative values)
            for counts in self.product_element_counts:
                row.append(-counts.get(element, 0))
            matrix.append(row)

        return matrix

    def find_positive_integer_solution(self):
        # Convert the coefficient matrix to a NumPy array
        A = np.array(self.coefficient_matrix)

        # Create a linear programming problem
        problem = LpProblem("Chemical_Equation_Balancer", sense=1)  # Minimize
        
        # Create variables for coefficients (all must be positive integers)
        num_variables = len(self.reactants) + len(self.products)
        coeffs = LpVariable.dicts("coeff", range(num_variables), lowBound=1, cat=LpInteger)
        
        # Add constraints for each element's balance equation
        for row in A:
            problem += lpSum(row[i] * coeffs[i] for i in range(num_variables)) == 0
            
        # Solve the problem
        problem.solve()
        
        if LpStatus[problem.status] == 'Optimal':
            return [int(coeffs[i].varValue) for i in range(num_variables)]
        
        return None

# Example usage
equation = Equation("AlO3H3 + H2SO4 = Al2S3O12 + H2O")
solution_coeffs = equation.find_positive_integer_solution()

if solution_coeffs:
    print("Positive integer solution coefficients:", solution_coeffs)
else:
    print("No valid solution found")

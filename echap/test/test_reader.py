from alinea.echap.interception_leaf import pesticide_applications

data = """date,dose, product_name
2000-10-02 01:00:00, 1.5, Opus
2000-10-01 05:00:00, 2, Banko 500
2000-10-01 08:00:00, 1.2, Opus
"""

calendar = pesticide_applications(data)
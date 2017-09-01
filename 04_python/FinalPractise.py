import sys
import PairwiseAlign
from flask import Flask, render_template, request, url_for

app = Flask(__name__)


@app.route('/')
def script():
    return render_template('PairwiseAlign_Template.html')


@app.route('/results', methods=['GET', 'POST'])
def results():
	input_seq1 = request.form['input_seq1']
	input_seq2 = request.form['input_seq2']
	M = PairwiseAlign.MGenerator(input_seq1, input_seq2)
	Align = PairwiseAlign.dinamicAnalizer(M)
	out_text1 = (Align[0][0])
	out_text2 = (Align[0][1])
	out_text3 = (Align[1])
	return render_template('AlignResults.html', output_message1=out_text1,
						   output_message2=out_text2,output_message3=out_text3)


app.run(debug=True)
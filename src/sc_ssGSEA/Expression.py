from abc import ABC, abstractmethod
import inspect
import math
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import spmatrix, csr_matrix
from typing import List, Optional, Type, Union
#import h5py
import scanpy as sc
import subprocess
import sys

#from .SeuratObject import SeuratObjectRDS
from .sc_ssGSEA import clr


class Expression(ABC):
	_filepath: str
	_group_name: str
	_chip_path: Optional[str]
	_metacells: Optional[pd.DataFrame]
	_gene_names: Optional[List[str]] 
	_cell_names: Optional[List[str]]
	_group_labels: Optional[pd.Series]

	def __init__(
		self,
		filepath: str,
		group_column_name: str,
		chip_file_path: Optional[str] = None
	) -> None:
		"""
		"""
		self._filepath = filepath
		self._group_name = group_column_name
		self._chip_path = chip_file_path

		self._metacells = None
		self._gene_names = None
		self._cell_names = None
		self._group_labels = None

	@abstractmethod
	def load(self) -> None:
		pass

	def _normalize_sparse_matrix(
		self,
		sparse_mat: csr_matrix,
		block_size = 5000
	) -> csr_matrix:
		"""
		Helper function for `load()`, wraps the `clr()` function
		"""
		new_data = []
		new_indices_col_indx = []
		new_indptr_row_indx = [0]

		block_starts = np.arange(0, sparse_mat.shape[0], block_size)

		print("Normalizing genes in blocks.")

		for block_start in block_starts:
			print(f"Normalizing genes {block_start} - {block_start + block_size} out of {sparse_mat.shape[0]}")

			dense_data = sparse_mat[block_start:block_start+block_size,:].toarray()

			for full_i in range(block_start, min(block_start+block_size, sparse_mat.shape[0])):

				block_i = full_i - block_start

				dense_row = dense_data[block_i,:]

				if sum(dense_row) == 0:
					new_indptr_row_indx.append(new_indptr_row_indx[-1])
					continue

				norm_dense_row = clr(dense_row)

				n_in_row = 0

				for j in range(len(norm_dense_row)):

					if norm_dense_row[j] != 0:
						new_data.append(norm_dense_row[j])
						new_indices_col_indx.append(j)
						n_in_row += 1

				new_indptr_row_indx.append(new_indptr_row_indx[-1] + n_in_row)

		#print(f"CSR new_data: {len(new_data)}")
		#print(f"CSR new_indices: {len(new_indices_col_indx)}")
		#print(f"CSR new_indptr: {len(new_indptr_row_indx)}")
		#print(f"CSR indptr last: {new_indptr_row_indx[-1]}")

		#new_indptr_row_indx.append(len(new_data)-1)

		return csr_matrix((new_data, new_indices_col_indx, new_indptr_row_indx))

		#return sparse_mat

	def _get_metacells(
		self,
		sparse_mat: csr_matrix
	) -> pd.DataFrame:
		"""
		Helper function for `load()`, replaces RDS-specific `get_metacells()`
		from `sc_ssGSEA.py`
		"""
		try:
			assert self._group_labels is not None
			assert self._gene_names is not None
			assert self._cell_names is not None
		except AssertionError:
			raise RuntimeError((
				"Missing one or more required fields _group_labels, _gene_names, _cell_names. "
				"If using a custom parser, ensure that all these fields are populated prior "
				"to calling `_get_metacells()`."
			))

		levels = self._group_labels.unique()

		mc_d = {}

		for level in levels:
			colnames_with_level = self._group_labels[self._group_labels == level].index.tolist()

			level_inds = [
				self._group_labels.index.tolist().index(colname)
				for colname in colnames_with_level
			]

			mc_d[f"{level}"] = sparse_mat[:,level_inds].T.mean(axis = 0).tolist()[0]

		mc_df = pd.DataFrame(mc_d, index = self._gene_names)

		return mc_df

	def _convert_gene_names(self):
		"""
		"""
		if self._gene_names is None:
			raise RuntimeError((
				"Parser cannot call _convert_gene_names() until "
				"gene names have been parsed." 
			))

		if self._chip_path is None:
			raise RuntimeError((
				"Parser cannot call _convert_gene_names() without "
				"providing a CHIP file."
			))

		chip_df = pd.read_csv(self._chip_path, sep = '\t')

		chip_map = {
			probe_id: symbol
			for probe_id, symbol in zip(chip_df["Probe Set ID"], chip_df["Gene Symbol"])
		}

		new_gene_names = []

		for gene in self._gene_names:
			try:
				new_gene_names.append(chip_map[gene])
			except KeyError:
				new_gene_names.append(gene)

		self._gene_names = new_gene_names

	@classmethod
	def get_expression_object(
		cls,
		filepath: str,
		group_column_name: str,
		chip_file_path: Optional[str] = None,
		custom_parser: Optional[Type["Expression"]] = None
	) -> "Expression":
		"""
		"""
		if custom_parser is not None:
			return custom_parser(
				filepath, 
				group_column_name,
				chip_file_path = chip_file_path
			)

		suffix = filepath.split('/')[-1].split('.')[-1]

		if suffix == "rds":
			return SeuratRDS_Expression(
				filepath, 
				group_column_name,
				chip_file_path = chip_file_path
			)
		elif suffix == "h5seurat":
			return H5Seurat_Expression(
				filepath, 
				group_column_name,
				chip_file_path = chip_file_path
			)
		elif suffix == "h5ad":
			return H5AD_Expression(
				filepath, 
				group_column_name,
				chip_file_path = chip_file_path
			)
		else:
			raise ValueError(f"No registered parser for suffix '{suffix}'.")

	@property
	def group_labels(self) -> pd.Series:
		"""
		"""
		if self._group_labels is None:
			raise RuntimeError(
				"Can't access group labels, have you called load()?"
			)
		else:
			return self._group_labels

	@property
	def metacells(self) -> pd.DataFrame:
		"""
		"""
		if self._metacells is None:
			raise RuntimeError(
				"Can't access metacells, have you called load()?"
			)
		else:
			return self._metacells


class SeuratRDS_Expression(Expression):

	def _call_seurat(self) -> None:
		module = inspect.getmodule(self)
		assert module is not None
		assert module.__file__ is not None
		script_path = '/'.join(module.__file__.split('/')[:-1]) + "/load_seurat.R"

		try:
			##TODO: FIX PATH AFTER SETTING PACKAGE SCRIPTS

			output = subprocess.check_output(
				["Rscript", "--vanilla", script_path, self._filepath],
				stderr = subprocess.STDOUT
			).decode()
		except subprocess.CalledProcessError as exc:
			output = exc.output.decode()

			if "there is no package" in output:
				raise RuntimeError((
					f"ERROR: Seurat was not found. To input an RDS file, please install "
					f"Seurat in this environment.\nSee instructions at "
					f"https://satijalab.org/seurat/articles/install.html\n"
					f"If you believe this message was generated in error, please "
					f"open an issue thread at https://github.com/genepattern/scGSEA/issues"
				))

			else:
				raise RuntimeError((
					f"ERROR: Unknown error when calling Seurat. Please make sure you "
					f"have Seurat installed.\nYou may report bugs at "
					f"https://github.com/genepattern/scGSEA/issues"
					f"\n Error text: {output}"
				))

		except FileNotFoundError:
			raise RuntimeError((
				f"ERROR: R not found. To input an RDS file, please install R and "
				f"the Seurat package.\n\tR installation: https://cran.r-project.org/mirrors.html\n\t" 
				f"Seurat installation: https://satijalab.org/seurat/articles/install.html\n"
				f"If you believe this message was generated in error, please "
				f"open an issue thread at https://github.com/genepattern/scGSEA/issues"
			))


	def load(self) -> None:
		self._call_seurat()

		with open("/tmp/genes.tsv", 'r') as f:
			self._gene_names = [
				line.strip('\n') for line in f.readlines()
			]

		with open("/tmp/barcodes.tsv", 'r') as f:
			self._cell_names = [
				line.strip('\n') for line in f.readlines()
			]

		sparse_mat = self._normalize_sparse_matrix(
			mmread("/tmp/expr.mtx").tocsr()
		)

		metadata_df = pd.read_csv(
			"/tmp/metadata.tsv",
			sep = '\t'
		)

		self._group_labels = metadata_df[self._group_name]

		self._metacells = self._get_metacells(sparse_mat)


class H5Seurat_Expression(Expression):
	def load(self) -> None:
		"""
		"""
		h5_file = h5py.File(self._filepath)
		active_assay_name = h5_file.attrs["active.assay"][0]

		## Extract expression and convert to sparse matrix + normalize

		expr = h5_file["assays"][active_assay_name]["counts"]

		sparse_mat = csr_matrix(
			(expr["data"], expr["indices"], expr["indptr"])
		).T

		sparse_mat = self._normalize_sparse_matrix(sparse_mat)

		## Extract row and column names

		self._gene_names = [
			gene.decode("utf-8") 
			for gene in h5_file["assays"][active_assay_name]["meta.features"]["_index"]
		]

		if self._chip_path is not None:
			self._convert_gene_names()

		self._cell_names = [
			barcode.decode("utf-8")
			for barcode in h5_file["cell.names"][:]
		]

		## Extract metadata and group labels

		meta_data = h5_file["meta.data"]

		self._group_labels = pd.Series(
			meta_data[self._group_name]["values"],
			index = meta_data["_index"]
		)

		self._metacells = self._get_metacells(sparse_mat)


class H5AD_Expression(Expression):
	def _probably_all_ints(
		self,
		vals: List[float]
	) -> bool:
		return all([
			math.isclose(x, round(x))
			for x in vals
		])

	def _mat_is_valid(
		self,
		X: np.ndarray | spmatrix | None
	) -> bool:
		"""
		"""
		if X is None:
			return False

		if isinstance(X, np.ndarray):
			first_row: List[float] = X[0,:].tolist()

		else: ## assume sparse
			first_row: List[float] = X[0,:].toarray()[0,:].tolist() #type: ignore

		return self._probably_all_ints(first_row)

	def _X_to_csr(
		self,
		X: np.ndarray | spmatrix | None
	) -> csr_matrix:
		"""
		"""
		assert X is not None

		if isinstance(X, np.ndarray):
			return csr_matrix(X).T

		else:
			return X.tocsr().T


	def load(self) -> None:
		"""
		"""

		## Extract expression from adata.raw.X (because must be counts)

		adata = sc.read_h5ad(self._filepath)
		raw_adata = adata.raw.to_adata()

		if not self._mat_is_valid(raw_adata.X):
			raise ValueError((
				"No counts found in input AnnData object. Object must have "
				".raw` populated, and data retrieved by calling "
				"`adata.raw.to_adata().X` must consist of values within a "
				"negligble rounding distance from zero (tolerance approximately "
				"1e-9). If your upstream analysis saves an explicit copy "
				"to `adata.raw`, that call should occur before any normalization "
				"or other transformation."
			))

		sparse_mat = self._X_to_csr(raw_adata.X)

		sparse_mat = self._normalize_sparse_matrix(sparse_mat)

		## Get gene names

		self._gene_names = raw_adata.var_names.to_list()

		## Get cell names

		self._cell_names = raw_adata.obs_names.to_list()

		## Get metadata

		self._group_labels = adata.obs.loc[:,self._group_name]

		if self._group_labels.shape[0] != len(self._cell_names):
			print((
				"Number of grouping labels from `adata.obs` and barcodes in "
				"`adata.raw.obs_names` do not match. Only using cells "
				"that have grouping information present in `adata.obs`"
			))

			shared_barcodes = list(set(self._group_labels.index).intersection(
				self._cell_names
			))

			self._group_labels = self._group_labels[shared_barcodes]
			self._cell_names = shared_barcodes

		## Get metacells

		self._metacells = self._get_metacells(sparse_mat)








